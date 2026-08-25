const RankT = Int64     # position of a simplex in its dimension's filtration order
const SlotT = Int64
const CombT = Int64     # combinatorial index of a simplex
const ValueT = Float32

# ------------------------------------------------------------
# FastMap: RankT key (must be >= 1) -> SlotT value
# Key -1 marks an empty slot
# Used by SetDict for R and Ri, by the active registry, and by FaceSlot when sparse
# ------------------------------------------------------------

mutable struct FastMap
    keys::Vector{RankT}
    vals::Vector{SlotT}
    mask::Int
    shift::Int      # 64 - log2(capacity): Fibonacci hashing takes the top bits
    count::Int
    limit::Int      
end

function FastMap(capacity::Int = 16)
    cap = 16
    while cap < capacity
        cap <<= 1
    end
    return FastMap(fill(RankT(-1), cap), zeros(SlotT, cap), cap - 1,
                   64 - trailing_zeros(cap), 0, (cap * 7) >> 3)
end

# Widen through Int64 first so this is correct for RankT = Int32 and Int64
@inline fm_hash(k::RankT, shift::Int) =
    Int((reinterpret(UInt64, Int64(k)) * 0x9E3779B97F4A7C15) >> shift) # Fibonacci hashing

# Get slot for key k
@inline function fm_get(m::FastMap, k::RankT)::SlotT
    i = fm_hash(k, m.shift)
    @inbounds while true
        kk = m.keys[i + 1]
        kk == k && return m.vals[i + 1]
        kk == -1 && return SlotT(0)
        i = (i + 1) & m.mask
    end
end

# Grow when too full
function fm_grow!(m::FastMap)
    oldk = m.keys
    oldv = m.vals
    cap = (m.mask + 1) << 1
    m.keys = fill(RankT(-1), cap)
    m.vals = zeros(SlotT, cap)
    m.mask = cap - 1
    m.shift = 64 - trailing_zeros(cap)
    m.limit = (cap * 7) >> 3
    @inbounds for t in eachindex(oldk)
        k = oldk[t]
        k == -1 && continue
        i = fm_hash(k, m.shift)
        while m.keys[i + 1] != -1
            i = (i + 1) & m.mask
        end
        m.keys[i + 1] = k
        m.vals[i + 1] = oldv[t]
    end
    return nothing
end

# Set a value
@inline function fm_set!(m::FastMap, k::RankT, v::SlotT)
    i = fm_hash(k, m.shift)
    @inbounds while true
        kk = m.keys[i + 1]
        if kk == k
            m.vals[i + 1] = v
            return nothing
        elseif kk == -1
            m.keys[i + 1] = k
            m.vals[i + 1] = v
            m.count += 1
            m.count > m.limit && fm_grow!(m)
            return nothing
        end
        i = (i + 1) & m.mask
    end
end

# returns k's value if present else 0, plus the index k belongs at
@inline function fm_find_or_hole(m::FastMap, k::RankT)
    i = fm_hash(k, m.shift)
    @inbounds while true
        kk = m.keys[i + 1]
        kk == k && return (m.vals[i + 1], i)
        kk == -1 && return (SlotT(0), i)
        i = (i + 1) & m.mask
    end
end

# Insert at a hole from fm_find_or_hole. Writes before growing, so the hole stays valid
@inline function fm_insert_at!(m::FastMap, hole::Int, k::RankT, v::SlotT)
    @inbounds m.keys[hole + 1] = k
    @inbounds m.vals[hole + 1] = v
    m.count += 1
    m.count > m.limit && fm_grow!(m)
    return nothing
end

# Backward-shift deletion: keeps the table tombstone-free so lookups never degrade after a long run of births and deaths
function fm_delete!(m::FastMap, k::RankT)
    i = fm_hash(k, m.shift)
    @inbounds while true
        kk = m.keys[i + 1]
        kk == -1 && return false
        kk == k && break
        i = (i + 1) & m.mask
    end
    mask = m.mask
    shift = m.shift
    j = i
    @inbounds while true
        j = (j + 1) & mask
        kj = m.keys[j + 1]
        if kj == -1
            break
        end
        h = fm_hash(kj, shift)
        if ((j - h) & mask) >= ((j - i) & mask)
            m.keys[i + 1] = kj
            m.vals[i + 1] = m.vals[j + 1]
            i = j
        end
    end
    @inbounds m.keys[i + 1] = -1
    m.count -= 1
    return true
end

# ------------------------------------------------------------
# SetDict: RankT key -> sorted Vector{RankT}
# ------------------------------------------------------------

# FastMap: Key --> slot 
# sets[slot] is corresponding value
# free stores slots not in use, if nonempty, pop a slot to store a new value (if empty grow sets)

# sing is vect, sing[i]=0 means sets[i] is not a singleton. Otherwise sing[i]=x means sets[i]={x}
# Benefit is sing is generally small enough to store in cache
mutable struct SetDict
    map::FastMap
    sets::Vector{Vector{RankT}}
    free::Vector{SlotT}
    sing::Vector{RankT}
    side_vpc::Int
    sdiam::Vector{ValueT}
    svert::Vector{Int32}
end

SetDict() = SetDict(FastMap(64), Vector{RankT}[], SlotT[], RankT[],
                    0, ValueT[], Int32[])

# update sing after sets[slot] has changed
@inline function sd_touch_shadow!(sd::SetDict, slot::SlotT)
    v = @inbounds sd.sets[slot]
    @inbounds sd.sing[slot] = length(v) == 1 ? v[1] : RankT(0)
    return nothing
end

@inline sd_slot(sd::SetDict, k::RankT)::SlotT = fm_get(sd.map, k) # key -> slot
@inline sd_set(sd::SetDict, slot::SlotT) = @inbounds sd.sets[slot] # slot -> sets[slot]
@inline sd_haskey(sd::SetDict, k::RankT)::Bool = fm_get(sd.map, k) != 0 # Check if key is a key

# Creates (or clears) the entry for k and returns its slot
function sd_create_slot!(sd::SetDict, k::RankT)::SlotT
    slot, hole = fm_find_or_hole(sd.map, k)
    if slot != 0 # defensive, k already a key, clear its slot and return
        empty!(@inbounds sd.sets[slot])
        @inbounds sd.sing[slot] = RankT(0)
        return slot
    end
    if isempty(sd.free)
        push!(sd.sets, RankT[])
        push!(sd.sing, RankT(0))
        if sd.side_vpc != 0
            push!(sd.sdiam, zero(ValueT))
            resize!(sd.svert, length(sd.sets) * sd.side_vpc)
        end
        slot = SlotT(length(sd.sets))
    else
        slot = pop!(sd.free)
        empty!(@inbounds sd.sets[slot])
    end
    @inbounds sd.sing[slot] = RankT(0)
    fm_insert_at!(sd.map, hole, k, slot)
    return slot
end

# Creates the entry k = {v}
function sd_singleton!(sd::SetDict, k::RankT, v::RankT)::SlotT
    slot = sd_create_slot!(sd, k) # create the vector for k, keeping its slot
    push!(@inbounds(sd.sets[slot]), v) # add v
    @inbounds sd.sing[slot] = v # update sing
    return slot
end

# record a key's own diameter and vertex list against its pool slot
@inline function sd_set_side!(sd::SetDict, slot::SlotT, dia::ValueT,
                              verts::Vector{Int32}, base::Int)
    sd.side_vpc == 0 && return nothing
    @inbounds sd.sdiam[slot] = dia
    dst = (Int(slot) - 1) * sd.side_vpc
    @inbounds for i in 1:sd.side_vpc
        sd.svert[dst + i] = verts[base + i]
    end
    return nothing
end

# rank -> its stored diameter / vertex-list base for live keys only
@inline sd_side_diam(sd::SetDict, k::RankT)::ValueT =
    @inbounds sd.sdiam[fm_get(sd.map, k)]
@inline sd_side_base(sd::SetDict, k::RankT)::Int =
    (Int(fm_get(sd.map, k)) - 1) * sd.side_vpc

# remove y from sets[slot]
@inline function sd_remove_elem!(sd::SetDict, slot::SlotT, y::RankT)
    sset_remove!(@inbounds(sd.sets[slot]), y)
    sd_touch_shadow!(sd, slot)
    return nothing
end

# sets[slot] <- sets[slot] XOR src, by merging into `scratch` and then
# swapping the two vectors rather than copying the result back
@inline function sd_xor_into!(sd::SetDict, slot::SlotT, src::Vector{RankT},
                              scratch::Vector{RankT})::Vector{RankT}
    dst = @inbounds sd.sets[slot]
    isempty(src) && return scratch
    if isempty(dst)
        resize!(dst, length(src))
        copyto!(dst, src)
        sd_touch_shadow!(sd, slot)
        return scratch
    end
    sset_xor2!(scratch, dst, src)
    @inbounds sd.sets[slot] = scratch
    sd_touch_shadow!(sd, slot)
    return dst
end

# delete k from dict
function sd_delete!(sd::SetDict, k::RankT)
    slot = fm_get(sd.map, k)
    slot == 0 && return nothing
    fm_delete!(sd.map, k)
    empty!(@inbounds sd.sets[slot])
    @inbounds sd.sing[slot] = RankT(0)
    push!(sd.free, slot)
    return nothing
end


# fast check done during active enumeration to avoid the entire xor if it would become zero (even if a face is active)
# observation is that the majority of the time this happens, it is something like {a}xor{a} with singletons
# this check determines if you would just get cancellations via singletons
# it is fast since sing is generally small and fast to access
# xslot = slot of x, the simplex active enumeration is run on
# beneficial to run before coface_ready to avoid the expensive check that x is the youngest face when possible too

@inline function seeded_is_noop(sd::SetDict, xslot::SlotT, seeds::Vector{SlotT})::Bool
    cnt = 0
    f = RankT(0)
    if xslot != 0
        @inbounds f = sd.sing[xslot] 
        f == 0 && return false # not a singleton
        cnt = 1
    end
    @inbounds for t in eachindex(seeds)
        v = sd.sing[seeds[t]]
        v == 0 && return false
        if cnt == 0
            f = v
        elseif v != f
            return false
        end
        cnt += 1
    end
    return iseven(cnt)
end

# Iterates the live keys 
function sd_keys!(out::Vector{RankT}, sd::SetDict)
    empty!(out)
    ks = sd.map.keys
    @inbounds for i in eachindex(ks)
        k = ks[i]
        k != -1 && push!(out, k)
    end
    return out
end

# ------------------------------------------------------------
# Sorted-vector set operations (Vector{Int32}, ascending)
# ------------------------------------------------------------
# remove x from v
@inline function sset_remove!(v::Vector{RankT}, x::RankT)
    i = searchsortedfirst(v, x)
    @inbounds if i <= length(v) && v[i] == x
        deleteat!(v, i)
    end
    return v
end
# toggle x in v
@inline function sset_toggle!(v::Vector{RankT}, x::RankT)
    i = searchsortedfirst(v, x)
    @inbounds if i <= length(v) && v[i] == x
        deleteat!(v, i)
    else
        insert!(v, i, x)
    end
    return v
end

# dst = dst XOR src using sorted sets
# add the smaller element and advance that set, on a tie advance both and add nothing
function sset_xor!(dst::Vector{RankT}, src::Vector{RankT}, buf::Vector{RankT})
    ni = length(dst)
    nj = length(src)
    nj == 0 && return dst
    if ni == 0
        resize!(dst, nj)
        copyto!(dst, src)
        return dst
    end
    length(buf) < ni + nj && resize!(buf, ni + nj)
    i = 1; j = 1; o = 0
    @inbounds while i <= ni && j <= nj
        a = dst[i]; b = src[j]
        if a < b
            o += 1; buf[o] = a; i += 1
        elseif a > b
            o += 1; buf[o] = b; j += 1
        else
            i += 1; j += 1
        end
    end
    @inbounds while i <= ni
        o += 1; buf[o] = dst[i]; i += 1
    end
    @inbounds while j <= nj
        o += 1; buf[o] = src[j]; j += 1
    end
    resize!(dst, o)
    @inbounds copyto!(dst, 1, buf, 1, o)
    return dst
end

# sorted set xor but output is a seperate vector so no buffer needs to be used
function sset_xor2!(out::Vector{RankT}, a::Vector{RankT}, b::Vector{RankT})
    na = length(a); nb = length(b)
    length(out) < na + nb && resize!(out, na + nb)
    i = 1; j = 1; o = 0
    @inbounds while i <= na && j <= nb
        x = a[i]; y = b[j]
        if x < y
            o += 1; out[o] = x; i += 1
        elseif x > y
            o += 1; out[o] = y; j += 1
        else
            i += 1; j += 1
        end
    end
    @inbounds while i <= na
        o += 1; out[o] = a[i]; i += 1
    end
    @inbounds while j <= nb
        o += 1; out[o] = b[j]; j += 1
    end
    resize!(out, o)
    return out
end

# ------------------------------------------------------------
# UnionFind for dimension 0
#
# R_0[v] is always the singleton {rep(v)} and Ri_0[r] is 
# r's connected component

# equivalent to union find where root is forced to be older 
# ------------------------------------------------------------

# parent[v] is v's parent
mutable struct UnionFind
    parent::Vector{Int32}
end

UnionFind(n::Int) = UnionFind(Int32[i for i in 1:n])

# find the root of x
@inline function uf_find(uf::UnionFind, x::Int32)::Int32
    p = uf.parent
    r = x
    @inbounds while p[r] != r # first pass finding root
        r = p[r]
    end
    @inbounds while p[x] != r # second pass updating path
        nx = p[x]
        p[x] = r
        x = nx
    end
    return r
end

# ------------------------------------------------------------
# FaceSlot: combinatorial index -> reduction-set pool slot (0 = not present)
# faces come as vertex sets, so this gets R[face] by going 
# comb index -> slot -> pool[slot] = R[face]

# also used by active enumeration, given a face gets 
# {set of all active n-simps contianin that face} via 
# comb index -> slot -> 
# ar.lr[slot]  = ranks of actives containing it
# ar.lw[slot]  = each one's extra vertex
# ar.ls[slot]  = each one's R-pool slot

# when dense gets too large, switches to sparse
# ------------------------------------------------------------

const DENSE_INDEX_CAP = 1 << 24

# ------------------------------------------------------------
# FaceSlot: combinatorial index -> pool slot (0 = not present)
# ------------------------------------------------------------

mutable struct FaceSlot
    dense::Vector{SlotT}
    sparse::FastMap
    isdense::Bool
end

function FaceSlot(universe::Int)
    if universe >= 0 && universe <= DENSE_INDEX_CAP
        return FaceSlot(zeros(SlotT, universe), FastMap(16), true)
    end
    return FaceSlot(SlotT[], FastMap(16), false)
end

@inline function fs_get(r::FaceSlot, k::CombT)::SlotT
    if r.isdense
        return @inbounds r.dense[k + 1] # offset since comb starts at 0
    end
    return fm_get(r.sparse, k)
end

@inline function fs_set!(r::FaceSlot, k::CombT, v::SlotT)
    if r.isdense
        @inbounds r.dense[k + 1] = v
    elseif v == 0 # 0 = delete
        fm_delete!(r.sparse, k)
    else
        fm_set!(r.sparse, k, v)
    end
    return nothing
end

# ------------------------------------------------------------
# Distance matrix 
# ------------------------------------------------------------

struct LowerDistances
    n::Int
    v::Vector{ValueT}
end

@inline lower_offset(i::Int, j::Int) = (((j - 1) * (j - 2)) >> 1) + i

Base.@propagate_inbounds @inline function Base.getindex(m::LowerDistances,
                                                        i::Int, j::Int)::ValueT
    i == j && return zero(ValueT)
    return i < j ? m.v[lower_offset(i, j)] : m.v[lower_offset(j, i)]
end

function LowerDistances(mat::AbstractMatrix)
    n = size(mat, 1)
    v = Vector{ValueT}(undef, (n * (n - 1)) >> 1)
    @inbounds for j in 2:n, i in 1:(j - 1)
        v[lower_offset(i, j)] = ValueT(mat[i, j])
    end
    return LowerDistances(n, v)
end

struct DistanceMatrix
    n::Int # num verts
    mat::LowerDistances
end

function DistanceMatrix(mat::AbstractMatrix)
    n = size(mat, 1)
    @assert size(mat, 2) == n "distance matrix must be square"
    return DistanceMatrix(n, LowerDistances(mat))
end

function DistanceMatrix(mat::LowerDistances)
    return DistanceMatrix(mat.n, mat)
end

# ------------------------------------------------------------
# Binomial table 
# ------------------------------------------------------------

struct BinomialTable
    data::Matrix{CombT}
end

function BinomialTable(nmax::Int, kmax::Int)
    data = zeros(CombT, kmax + 1, nmax + 1)
    @inbounds for n in 0:nmax
        data[1, n + 1] = 1
        for k in 1:min(n, kmax)
            if k == n
                data[k + 1, n + 1] = 1
            else
                data[k + 1, n + 1] = data[k, n] + data[k + 1, n]
            end
        end
    end
    return BinomialTable(data)
end

@inline function choose(binom::BinomialTable, n::Int, k::Int)::CombT
    if k > n || k < 0 || n < 0
        return CombT(0)
    end
    return @inbounds binom.data[k + 1, n + 1]
end

# Combinatorial index of an ascending vertex list of length m
@inline function comb_index(verts::Vector{Int32}, base::Int, m::Int, binom::BinomialTable)::CombT
    idx = CombT(0)
    @inbounds for k in 1:m
        idx += choose(binom, Int(verts[base + k]), k)
    end
    return idx
end

# Combinatorial index of the face obtained by omitting position omit
@inline function face_index_without(verts::Vector{Int32}, base::Int, m::Int, omit::Int,
                                    binom::BinomialTable)::CombT
    idx = CombT(0)
    @inbounds begin
        for k in 1:(omit - 1)
            idx += choose(binom, Int(verts[base + k]), k)
        end
        for k in (omit + 1):m
            idx += choose(binom, Int(verts[base + k]), k - 1)
        end
    end
    return idx
end

@inline function face_index_without1(verts::Vector{Int32}, m::Int, omit::Int,
                                     binom::BinomialTable)::CombT
    return face_index_without(verts, 0, m, omit, binom)
end

# ------------------------------------------------------------
# Bucket sort of a generated level into filtration order
# ------------------------------------------------------------

# find the bucket based on diameter d. lo = range min, inv_span = 1/(hi-lo), nb = num buckets
@inline function bucket_of(d::ValueT, lo::ValueT, inv_span::ValueT, nb::Int)::Int
    b = floor(Int, (d - lo) * inv_span * nb)
    b < 0 && return 0
    b >= nb && return nb - 1
    return b
end

# ------------------------------------------------------------
# EdgeList: the dimension-1 simplices in filtration order
# ------------------------------------------------------------

struct EdgeList
    count::Int
    diam::Vector{ValueT}        # diameters
    verts::Vector{Int32}        # vertices, flat, 2 per edge: edge r is [2r-1, 2r]
end

@inline ed_base(r::Integer) = (Int(r) - 1) * 2
@inline ed_diam(ed::EdgeList, r::Integer) = @inbounds ed.diam[Int(r)]
@inline ed_idx(ed::EdgeList, r::Integer, binom::BinomialTable) =
    comb_index(ed.verts, ed_base(r), 2, binom)

# ------------------------------------------------------------
# Sorting a level into filtration order 
# ------------------------------------------------------------

mutable struct SortScratch
    dbuf::Vector{ValueT}
    vbuf::Vector{Int32}
    cnt::Vector{RankT} # histogram of bucket counts
    sdbuf::Vector{ValueT}
    svbuf::Vector{Int32}
    skv::Vector{Tuple{ValueT,CombT,RankT}}
end

SortScratch() = SortScratch(ValueT[], Int32[], RankT[], ValueT[], Int32[],
                            Tuple{ValueT,CombT,RankT}[])

# sort a bucket, a = start index, b = end index
function sort_slice!(diam::Vector{ValueT}, verts::Vector{Int32}, vpc::Int,
                     a::Int, b::Int, binom::BinomialTable, sc::SortScratch, depth::Int)
    n = b - a + 1
    n <= 1 && return nothing # singleton

    if n <= 32 || depth >= 5 # small or depth >=5 directly sort
        sort_slice_small!(diam, verts, vpc, a, b, binom, sc, n > 32)
        return nothing
    end

    dmin = ValueT(Inf); dmax = ValueT(-Inf)
    @inbounds for i in a:b
        d = diam[i]
        d < dmin && (dmin = d)
        d > dmax && (dmax = d)
    end

    if dmax == dmin # recursive on diam offers nothing
        sort_slice_small!(diam, verts, vpc, a, b, binom, sc, true)
        return nothing
    end

    nb = n >> 2 # n/4 buckets
    nb < 2 && (nb = 2) # floor of two, unreachable but safe if n cap changes
    inv = ValueT(1) / (dmax - dmin) 
    length(sc.cnt) < nb + 1 && resize!(sc.cnt, nb + 1)
    cnt = sc.cnt
    @inbounds for t in 1:(nb + 1) # zero old counts
        cnt[t] = 0
    end
    @inbounds for i in a:b # make histogram
        cnt[bucket_of(diam[i], dmin, inv, nb) + 2] += RankT(1) 
    end
    @inbounds for t in 2:(nb + 1) #cnt[1] = 0, cnt[bk+1] = number of elements in buckets 0 â€¦ bk-1 = bucket bk's 0-based start offset
        cnt[t] += cnt[t - 1]
    end
    
    bounds = Vector{RankT}(undef, nb + 1)
    @inbounds for t in 1:(nb + 1)
        bounds[t] = cnt[t]
    end

    # resize diam and vertex buffs 
    length(sc.dbuf) < n && resize!(sc.dbuf, n)
    length(sc.vbuf) < n * vpc && resize!(sc.vbuf, n * vpc)

    # sort into buckets
    @inbounds for i in a:b
        d = diam[i]
        bk = bucket_of(d, dmin, inv, nb)
        p = cnt[bk + 1] + RankT(1)
        cnt[bk + 1] = p
        sc.dbuf[p] = d
        src = (i - 1) * vpc
        dst = (Int(p) - 1) * vpc
        for q in 1:vpc
            sc.vbuf[dst + q] = verts[src + q]
        end
    end

    # copy back
    @inbounds for t in 1:n
        diam[a + t - 1] = sc.dbuf[t]
        src = (t - 1) * vpc
        dst = (a + t - 2) * vpc
        for q in 1:vpc
            verts[dst + q] = sc.vbuf[src + q]
        end
    end

    # recursion on each new bucket
    @inbounds for bk in 1:nb
        sa = a + Int(bounds[bk])
        sb = a + Int(bounds[bk + 1]) - 1
        sb > sa && sort_slice!(diam, verts, vpc, sa, sb, binom, sc, depth + 1)
    end
    return nothing
end

# direct sort for n<=32, depth >=5, or if dim_dim==dim_max
# comparison sort used when n>32, else insertion sort
function sort_slice_small!(diam::Vector{ValueT}, verts::Vector{Int32}, vpc::Int,
                           a::Int, b::Int, binom::BinomialTable, sc::SortScratch,
                           use_comparison_sort::Bool)
    n = b - a + 1
    n <= 1 && return nothing
    length(sc.skv) < n && resize!(sc.skv, n)
    length(sc.sdbuf) < n && resize!(sc.sdbuf, n)
    length(sc.svbuf) < n * vpc && resize!(sc.svbuf, n * vpc)
    kv = sc.skv
    @inbounds for t in 1:n
        p = a + t - 1
        kv[t] = (diam[p], comb_index(verts, (p - 1) * vpc, vpc, binom), RankT(p)) # triples to sort
    end
    if use_comparison_sort
        sort!(view(kv, 1:n)) # built in comparison sort
    else
        # insertion sort
        @inbounds for i in 2:n
            x = kv[i]
            j = i - 1
            while j >= 1 && kv[j] > x
                kv[j + 1] = kv[j]
                j -= 1
            end
            kv[j + 1] = x
        end
    end
    dtmp = sc.sdbuf
    vtmp = sc.svbuf
    # order in buffer
    @inbounds for t in 1:n
        p = Int(kv[t][3])
        dtmp[t] = diam[p]
        src = (p - 1) * vpc
        dst = (t - 1) * vpc
        for q in 1:vpc
            vtmp[dst + q] = verts[src + q]
        end
    end
    # copy back
    @inbounds for t in 1:n
        diam[a + t - 1] = dtmp[t]
        src = (t - 1) * vpc
        dst = (a + t - 2) * vpc
        for q in 1:vpc
            verts[dst + q] = vtmp[src + q]
        end
    end
    return nothing
end

# lv_diam_v = diameters 
# lv_verts = vertices 
# vpc = vertices per simplex 
# starts = starts of buckets 
# nb = number of buckets
function finish_buckets!(lv_diam_v::Vector{ValueT}, lv_verts::Vector{Int32}, vpc::Int,
                         starts::Vector{RankT}, nb::Int, binom::BinomialTable,
                         sc::SortScratch)
    @inbounds for b in 1:nb
        a = Int(starts[b]) + 1
        e = Int(starts[b + 1])
        e > a && sort_slice!(lv_diam_v, lv_verts, vpc, a, e, binom, sc, 1)
    end
    return nothing
end

# ------------------------------------------------------------
# Dimension 1: edges, in filtration order
# first scan builds the diameter histogram, second puts edges in buckets
# ------------------------------------------------------------

function build_edges(dist::DistanceMatrix, threshold::ValueT, sc::SortScratch,
                     binom::BinomialTable)
    n = dist.n
    mat = dist.mat

    cap = (n * (n - 1)) >> 1
    nb = cap >> 2 # 4 edges per bucket
    nb < 1024 && (nb = 1024) # lower bound
    nb > (1 << 16) && (nb = 1 << 16) # upper bound for memory
    dmin = zero(ValueT)
    inv_span = (isfinite(threshold) && threshold > 0) ? ValueT(1) / threshold : ValueT(0)
    cnt = zeros(RankT, nb + 1)

    # pass 1: diameter histogram
    @inbounds for j in 1:(n - 1)
        for i in 0:(j - 1)
            d = mat[i + 1, j + 1]
            if d <= threshold
                cnt[bucket_of(d, dmin, inv_span, nb) + 2] += RankT(1)
            end
        end
    end

    @inbounds for b in 2:(nb + 1)
        cnt[b] += cnt[b - 1]
    end
    count = Int(cnt[nb + 1])      
    if count == 0
        return EdgeList(0, ValueT[], Int32[])
    end
    starts = copy(cnt)
    diam = Vector{ValueT}(undef, count)
    verts = Vector{Int32}(undef, 2 * count)

    # pass 2
    @inbounds for j in 1:(n - 1)
        for i in 0:(j - 1)
            d = mat[i + 1, j + 1]
            if d <= threshold
                b = bucket_of(d, dmin, inv_span, nb)
                p = cnt[b + 1] + RankT(1)
                cnt[b + 1] = p
                diam[p] = d
                verts[2 * Int(p) - 1] = Int32(i)
                verts[2 * Int(p)] = Int32(j)
            end
        end
    end

    finish_buckets!(diam, verts, 2, starts, nb, binom, sc)
    return EdgeList(count, diam, verts)
end

# ------------------------------------------------------------
# Dimensions >= 2: edge-driven streaming generation
# ------------------------------------------------------------

# g = matrix storing G_<e, the graph with all edges less than current edge
# need = number of vertices needed to form a d simplex when union with edge
# emit = function to process simplex
function colex_rec!(g::Matrix{UInt64}, words::Int, nv::Int, scratch::Matrix{UInt64},
                    need::Int, depth::Int, stack::Vector{Int32}, u::Int, v::Int,
                    merged::Vector{Int32}, vpc::Int, binom::BinomialTable,
                    single::Bool, akey::Vector{CombT}, avert::Vector{Int32},
                    dia::ValueT, emit::F) where {F}
    @inbounds for w in 1:words # number of words
        # this depth's candidate set, word by word 
        # candidate set = \cap of N(x) for all current verts, but only those verts less than the last vert added
        wb = scratch[w, depth]
        while wb != 0
            bpos = trailing_zeros(wb) # lowest bit position
            wb &= wb - 1 # clear that position
            x = (w - 1) * 64 + bpos # vertex idx
            x >= nv && break # safeguard against leftover higher bits, nothing ever sets those bits so never actually fired
            stack[depth] = Int32(x) # this depths vert
            if need == 1 
                # merge {u,v} with stack[1..depth] (descending) into ascending
                a = 1; b = depth; t = 0
                while a <= 2 && b >= 1
                    bv = (a == 1) ? Int32(u) : Int32(v)
                    if bv < stack[b]
                        t += 1; merged[t] = bv; a += 1
                    else
                        t += 1; merged[t] = stack[b]; b -= 1
                    end
                end
                while a <= 2
                    t += 1; merged[t] = (a == 1) ? Int32(u) : Int32(v); a += 1
                end
                while b >= 1
                    t += 1; merged[t] = stack[b]; b -= 1
                end
                if single # edge is alone at that length
                    emit(dia, merged, 0)
                else # multiple edges have same length, store comb index to process them with streams created by other edges
                    push!(akey, comb_index(merged, 0, vpc, binom))
                    for q in 1:vpc
                        push!(avert, merged[q])
                    end
                end
            else
                for t in 1:words
                    scratch[t, depth + 1] = scratch[t, depth] & g[t, x + 1] # intersection with x
                end
                cw = x >> 6 # keep only vertices below x
                for t in (cw + 2):words
                    scratch[t, depth + 1] = 0
                end
                rem = x & 63
                scratch[cw + 1, depth + 1] &= (rem == 0 ? UInt64(0) :
                                               ((UInt64(1) << rem) - 1))
                colex_rec!(g, words, nv, scratch, need - 1, depth + 1, stack, u, v,
                           merged, vpc, binom, single, akey, avert, dia, emit)
            end
        end
    end
    return nothing
end

function generate_stream(edges::EdgeList, nv::Int, d::Int, binom::BinomialTable,
                         emit::F) where {F}
    edges.count == 0 && return nothing
    words = cld(nv, 64)
    g = zeros(UInt64, words, nv)
    need = d - 1
    vpc = d + 1
    stack = zeros(Int32, need + 1)
    merged = zeros(Int32, vpc + 1)
    scratch = zeros(UInt64, words, need + 1)
    akey = CombT[]
    avert = Int32[]

    E = edges.count
    i = 1
    while i <= E
        r = @inbounds edges.diam[i]
        jend = i
        # need to walk in edge length tie groups
        @inbounds while jend <= E && edges.diam[jend] == r
            jend += 1
        end
        single = (jend == i + 1) # true if no ties, false otherwise
        empty!(akey); empty!(avert)
        for kk in i:(jend - 1)
            u = Int(@inbounds edges.verts[2 * kk - 1])
            v = Int(@inbounds edges.verts[2 * kk])
            anyb = UInt64(0)
            @inbounds for t in 1:words
                y = g[t, u + 1] & g[t, v + 1] # intersection of words
                scratch[t, 1] = y
                anyb |= y
            end
            if anyb != 0
                colex_rec!(g, words, nv, scratch, need, 1, stack, u, v, merged,
                           vpc, binom, single, akey, avert, r, emit)
            end
            # update g after edge has been processed
            @inbounds g[(v >> 6) + 1, u + 1] |= (UInt64(1) << (v & 63)) 
            @inbounds g[(u >> 6) + 1, v + 1] |= (UInt64(1) << (u & 63))
        end
        if !single && !isempty(akey) # sort and process all generated simps from ties
            ord = sortperm(akey)
            @inbounds for t in ord
                emit(r, avert, (t - 1) * vpc)
            end
        end
        i = jend
    end
    return nothing
end

#############################################################################
# RedZeD_VR algorithm 
#############################################################################

# ============================================================
# Active-simplex registry
#
# The keys of R_n are exactly the active n-simplices
# The registry keeps its vertex list
# for every (n-1)-face the list of active simplices containing it (with the extra vertex)
# ============================================================

mutable struct ActiveReg
    m::Int                          # vertices per top simplex (= n + 1)
    slot_of::FastMap                # active rank -> vertex-slot
    verts::Vector{Int32}            # flat, m per slot (VERTEX ids, not ranks)
    free::Vector{SlotT}
    nslots::Int
    face_slot::FaceSlot            # face comb index -> list slot
    lr::Vector{Vector{RankT}}       # per list slot: active ranks
    lw::Vector{Vector{Int32}}       # per list slot: the extra vertex
    ls::Vector{Vector{SlotT}}       # per list slot: that active's R-set pool slot
    lfree::Vector{SlotT}
end

function ActiveReg(m::Int, face_universe::Int)
    return ActiveReg(m, FastMap(1024), Int32[], SlotT[], 0,
                     FaceSlot(face_universe),
                     Vector{RankT}[], Vector{Int32}[], Vector{SlotT}[], SlotT[])
end

# allocate a new slot of free is empty, otherwise pop a slot
function ar_alloc_vslot!(ar::ActiveReg)::SlotT
    if isempty(ar.free)
        ar.nslots += 1
        resize!(ar.verts, ar.nslots * ar.m)
        return SlotT(ar.nslots)
    end
    return pop!(ar.free)
end

# register an active simplex that is active after active enumeration 
function register_active!(ar::ActiveReg, rank::RankT, rslot::SlotT, verts::Vector{Int32},
                          binom::BinomialTable)
    slot = ar_alloc_vslot!(ar)
    base = (Int(slot) - 1) * ar.m
    @inbounds for i in 1:ar.m
        ar.verts[base + i] = verts[i]
    end
    fm_set!(ar.slot_of, rank, slot)
    @inbounds for omit in 1:ar.m
        f = face_index_without(ar.verts, base, ar.m, omit, binom)
        ls = fs_get(ar.face_slot, f)
        if ls == 0
            if isempty(ar.lfree)
                push!(ar.lr, RankT[])
                push!(ar.lw, Int32[])
                push!(ar.ls, SlotT[])
                ls = SlotT(length(ar.lr))
            else
                ls = pop!(ar.lfree)
            end
            fs_set!(ar.face_slot, f, ls)
        end
        push!(ar.lr[ls], rank)
        push!(ar.lw[ls], ar.verts[base + omit])
        push!(ar.ls[ls], rslot)
    end
    return nothing
end

# doomed = set of ranks that just became inactive, remove them form active reg if present
function unregister_active!(ar::ActiveReg, doomed::Vector{RankT}, binom::BinomialTable)
    @inbounds for rank in doomed
        slot = fm_get(ar.slot_of, rank)
        slot == 0 && continue
        base = (Int(slot) - 1) * ar.m
        for omit in 1:ar.m
            f = face_index_without(ar.verts, base, ar.m, omit, binom)
            ls = fs_get(ar.face_slot, f)
            ls == 0 && continue
            lr = ar.lr[ls]
            lw = ar.lw[ls]
            lsl = ar.ls[ls]
            for k in eachindex(lr)
                if lr[k] == rank
                    lr[k] = lr[end]; pop!(lr)
                    lw[k] = lw[end]; pop!(lw)
                    lsl[k] = lsl[end]; pop!(lsl)
                    break
                end
            end
            if isempty(lr)
                fs_set!(ar.face_slot, f, SlotT(0))
                push!(ar.lfree, ls)
            end
        end
        fm_delete!(ar.slot_of, rank)
        push!(ar.free, slot)
    end
    return nothing
end

# Reusable Step-1 state
mutable struct SeedState
    seen::Vector{Bool} # seen verts as a vect
    seenbits::Vector{UInt64} # same info as a bitset to make step 2 subtraction fast
    verts::Vector{Int32} # set of candidate vertices
    seedlists::Vector{Vector{SlotT}} # R-set slots of the face-sharing actives
end

function SeedState(nv::Int, words::Int)
    return SeedState(fill(false, nv), zeros(UInt64, words), Int32[],
                     [SlotT[] for _ in 1:nv])
end

function reset_seeds!(s::SeedState)
    @inbounds for w in s.verts # only need to reset vertices actually found, why verts is stored
        s.seen[w + 1] = false
        empty!(s.seedlists[w + 1])
        s.seenbits[(Int(w) >> 6) + 1] = 0
    end
    empty!(s.verts)
    return nothing
end

# ============================================================
# Whole-run state
# ============================================================

mutable struct RState
    nv::Int # vertex count
    n::Int # homology degree 
    dist::DistanceMatrix # distance matrix
    binom::BinomialTable # binomial table 

    R::Vector{SetDict} # vect of R's as set dicts
    Ri::Vector{SetDict} # vect of Ri's as set dicts
    pairs::Vector{Vector{Tuple{ValueT,ValueT}}} # pairs found
    edges::EdgeList # sorted edges

    # face --> pool slot lookup
    csmap::FaceSlot
    cson::Bool
    csdim::Int

    uf::UnionFind # union find for dimension 0

    bar::Vector{RankT} # r \partial x (called bar because notation was origonally r\partial = \overline{\partial})
    buf::Vector{RankT} # buffer for xor
    touched::Vector{RankT} # snapshot of Ri[y]
    others::Vector{RankT} # bar \ {j} : entries whose Ri needs updating 
    removed::Vector{RankT} # keys removed from R
    topbar::Vector{RankT} # r \partial for n+1 dimensional simplices
    swapbuf::Vector{RankT} # buffer

    act::ActiveReg
    seeds::SeedState
    rbits::Matrix{UInt64} # bitset adjacency grown one edge at a time 
    rlo::Vector{Int32} # lowest bit for a vertex
    rhi::Vector{Int32} # highest bit for a vertex
    rptr::Int # pointer into edge list 
    words::Int # num words
    coface_buf::Vector{Int32} # buffer for n+1 simplex
    top_buf::Vector{Int32} # buffer for n simplex

    # optional simplex pairs: one entry per stored interval of pairs[k] in the same order
    sbirth::Vector{Vector{Int32}} # flat birth-simplex vertices, k per entry of pairs[k]
    sdeath::Vector{Vector{Int32}} # flat death-simplex vertices, k+1 per entry, all -1 = infinite
    want_spx::Bool # record simplex pairs, left empty when false
end

# ============================================================
# Reduction bookkeeping
# ============================================================

# get the diameter of a p dimensional simplex r 
# Ri used for diams since only ever called on alive elements
@inline function lvl_diam_of(st::RState, p::Int, r::RankT)::ValueT
    p == 1 && return ed_diam(st.edges, r)
    return sd_side_diam(st.Ri[p], r)
end

# get combinatorial index of a p dimensional simplex r 
# side base is offset into the list of vertices
@inline function lvl_idx_of(st::RState, p::Int, r::RankT)::CombT
    p == 1 && return ed_idx(st.edges, r, st.binom)
    sd = st.R[p]
    return comb_index(sd.svert, sd_side_base(sd, r), p + 1, st.binom)
end

# append the vertices of the p dimensional simplex r to out
@inline function append_verts!(out::Vector{Int32}, st::RState, p::Int, r::RankT)
    if p == 1
        b = ed_base(r)
        @inbounds push!(out, st.edges.verts[b + 1], st.edges.verts[b + 2])
    else
        sd = st.Ri[p]
        b = sd_side_base(sd, r)
        @inbounds for i in 1:(p + 1)
            push!(out, sd.svert[b + i])
        end
    end
    return nothing
end

# remove z from csmap (comb index --> slot)
# apply sometimes false since there isn't always a csmap (i.e. H1 or during AE)
@inline function cs_clear!(st::RState, z::RankT, apply::Bool)
    apply || return nothing
    fs_set!(st.csmap, lvl_idx_of(st, st.csdim, z), SlotT(0))
    return nothing
end

# special case for if r \partial = {y}, remove y from each R[x] st y in R[y]
function remove_maximal!(st::RState, Rd::SetDict, Rid::SetDict, y::RankT, apply_cs::Bool)
    removed = st.removed
    empty!(removed)
    uslot = sd_slot(Rid, y)
    touched = st.touched
    empty!(touched)
    append!(touched, sd_set(Rid, uslot)) # Ri[y]
    sd_delete!(Rid, y) # y dies 
    @inbounds for z in touched
        zslot = sd_slot(Rd, z) # never 0: z in Ri[y] means y in R[z], so z is a key of R
        sd_remove_elem!(Rd, zslot, y)
        if isempty(sd_set(Rd, zslot)) # z now inactive, remove
            cs_clear!(st, z, apply_cs)
            sd_delete!(Rd, z)
            push!(removed, z)
        end
    end
    return nothing
end

# general case, bar = r partial size geq 2, j = youngest element
function replace_pivot!(st::RState, Rd::SetDict, Rid::SetDict, bar::Vector{RankT}, j::RankT, apply_cs::Bool)
    removed = st.removed
    empty!(removed)
    uslot = sd_slot(Rid, j) # never 0: j = bar[end], same as remove_maximal!

    touched = st.touched
    empty!(touched)
    append!(touched, sd_set(Rid, uslot)) # Ri[j]

    others = st.others # others = bar \ {j}
    empty!(others)
    @inbounds for y in bar
        y == j && continue
        push!(others, y)
    end

    @inbounds for z in touched
        zslot = sd_slot(Rd, z) # never 0: z in Ri[j] means j in R[z], so z is a key of R
        st.swapbuf = sd_xor_into!(Rd, zslot, bar, st.swapbuf)
        if isempty(sd_set(Rd, zslot))
            cs_clear!(st, z, apply_cs)
            sd_delete!(Rd, z)
            push!(removed, z)
        end
    end
    buf = st.buf

    K = length(touched)
    @inbounds for y in others
        # never 0: y in bar means y is an element of some R set, so a key of Ri 
        sy = sd_set(Rid, sd_slot(Rid, y))
        if K * 128 < length(sy) # toggle is faster when touched is significalty smaller than sy
            for z in touched
                sset_toggle!(sy, z)
            end
        else
            sset_xor!(sy, touched, buf)
        end
        isempty(sy) && sd_delete!(Rid, y)
    end

    sd_delete!(Rid, j)
    return nothing
end

# handles r partial x \neq 0 for x of dimension d
# sdiam = simplex diam, used to store pair
function resolve_bar!(st::RState, d::Int, bar::Vector{RankT}, sdiam::ValueT)
    Rd = st.R[d - 1]
    Rid = st.Ri[d - 1]
    apply_cs = st.cson && (d - 1 == st.csdim)
    j = @inbounds bar[end]
    jdiam = lvl_diam_of(st, d - 1, j)

    # j is the birth simplex of the interval, recorded before its pool slot is freed
    # death simplex is the current d dimensional simplex 
    # left in coface_buf by coface_ready! in the top dimension and in top_buf by the generator below it
    if st.want_spx && sdiam != jdiam
        append_verts!(st.sbirth[d], st, d - 1, j)
        dv = (d == st.n + 1) ? st.coface_buf : st.top_buf
        @inbounds for k in 1:(d + 1)
            push!(st.sdeath[d], dv[k])
        end
    end

    if length(bar) == 1
        remove_maximal!(st, Rd, Rid, j, apply_cs)
    else
        replace_pivot!(st, Rd, Rid, bar, j, apply_cs)
    end

    if d == st.n + 1 && !isempty(st.removed) # top dimension, need to unregister no longer active 
        unregister_active!(st.act, st.removed, st.binom)
    end

    sdiam != jdiam && push!(st.pairs[d], (jdiam, sdiam)) # store pairs of nonzero length
    return nothing
end

# ============================================================
# Active enumeration
# ============================================================

# is b ~ a
@inline function rbit(st::RState, a::Int, b::Int)::Bool
    return @inbounds (st.rbits[(b >> 6) + 1, a + 1] >> (b & 63)) & 1 != 0
end

# advance rbits to new radius r
function advance_rbits!(st::RState, r::ValueT)
    lv = st.edges
    dm = lv.diam
    vt = lv.verts
    p = st.rptr
    cnt = lv.count
    rb = st.rbits
    rlo = st.rlo
    rhi = st.rhi
    @inbounds while p <= cnt && dm[p] <= r
        i = Int(vt[2 * p - 1])
        j = Int(vt[2 * p])
        wj = (j >> 6) + 1
        wi = (i >> 6) + 1
        rb[wj, i + 1] |= (UInt64(1) << (j & 63))
        rb[wi, j + 1] |= (UInt64(1) << (i & 63))
        wj32 = Int32(wj); wi32 = Int32(wi)
        rlo[i + 1] > wj32 && (rlo[i + 1] = wj32)
        rhi[i + 1] < wj32 && (rhi[i + 1] = wj32)
        rlo[j + 1] > wi32 && (rlo[j + 1] = wi32)
        rhi[j + 1] < wi32 && (rhi[j + 1] = wi32)
        p += 1
    end
    st.rptr = p
    return nothing
end

# add a vertex to a sortecd list of vertices
@inline function merge_vertex!(out::Vector{Int32}, xverts::Vector{Int32}, m::Int, w::Int32)::Int
    i = 1; j = 1
    @inbounds while i <= m && xverts[i] < w
        out[j] = xverts[i]; i += 1; j += 1
    end
    @inbounds out[j] = w
    wpos = j
    j += 1
    @inbounds while i <= m
        out[j] = xverts[i]; i += 1; j += 1
    end
    return wpos
end

# diameter of the face given by omitting the vertex omit 
@inline function face_diameter(verts::Vector{Int32}, m::Int, omit::Int, mat::LowerDistances)::ValueT
    diam = zero(ValueT)
    @inbounds for a in 1:m
        a == omit && continue
        va = Int(verts[a]) + 1
        for b in (a + 1):m
            b == omit && continue
            d = mat[va, Int(verts[b]) + 1]
            d > diam && (diam = d)
        end
    end
    return diam
end

# Is xverts U {w} ready to appear, i.e. is xverts its youngest face?
function coface_ready!(candverts::Vector{Int32}, xverts::Vector{Int32}, m::Int, w::Int32,
                       xidx::CombT, xdiam::ValueT, mat::LowerDistances, binom::BinomialTable)::Bool
    total = m + 1
    wpos = merge_vertex!(candverts, xverts, m, w)
    @inbounds for omit in 1:total
        omit == wpos && continue
        faceidx = face_index_without1(candverts, total, omit, binom)
        # xverts U {w} always has diam = diam xverts, so if faceidx<idx, face already born
        faceidx < xidx && continue 
        facedia = face_diameter(candverts, total, omit, mat)
        facedia < xdiam || return false # if the faceidx > idx, needs to have smaller diam
    end
    return true
end

# Step 1: candidate vertices w such that xverts U {w} is a valid coface containing  
# another active simplex
function collect_candidates!(st::RState, xverts::Vector{Int32})
    s = st.seeds
    reset_seeds!(s)
    ar = st.act
    m = ar.m
    binom = st.binom
    @inbounds for omit in 1:m # faces of x 
        f = face_index_without1(xverts, m, omit, binom)
        ls = fs_get(ar.face_slot, f)
        ls == 0 && continue
        xextra = Int(xverts[omit])
        lw = ar.lw[ls] # extra verts
        lsl = ar.ls[ls] # pool slots
        for t in eachindex(lw)
            w = lw[t]
            if rbit(st, xextra, Int(w)) # d(xextra,w)<=r, store 
                wi = Int(w) + 1
                if !s.seen[wi]
                    s.seen[wi] = true
                    push!(s.verts, w)
                    s.seenbits[(Int(w) >> 6) + 1] |= (UInt64(1) << (Int(w) & 63))
                end
                push!(s.seedlists[wi], lsl[t]) # add (x \ omit U {w}): by end of fcn stores all other active faces of x U w
            end
        end
    end
    return nothing
end

# Step 2: does another coface exist not yet processed
function has_leftover_coface!(st::RState, xverts::Vector{Int32}, m::Int, xidx::CombT,
                              r::ValueT)::Bool
    rb = st.rbits
    nv = st.nv
    sb = st.seeds.seenbits
    mat = st.dist.mat
    binom = st.binom
    cb = st.coface_buf
    # only words nonzero for all verts can survive intersection
    wlo = Int(0); whi = Int(typemax(Int32))
    @inbounds for i in 1:m
        v = Int(xverts[i]) + 1
        lo = Int(st.rlo[v]); hi = Int(st.rhi[v])
        lo > wlo && (wlo = lo)
        hi < whi && (whi = hi)
    end

    @inbounds for w in wlo:whi
        acc = rb[w, Int(xverts[1]) + 1]
        acc == 0 && continue
        for i in 2:m # repeated intersection
            acc &= rb[w, Int(xverts[i]) + 1]
            acc == 0 && break
        end
        acc == 0 && continue
        acc &= ~sb[w] # intersect with compliment of seen 
        while acc != 0
            bpos = trailing_zeros(acc)
            acc &= acc - 1
            cand = (w - 1) * 64 + bpos
            cand >= nv && break
            if coface_ready!(cb, xverts, m, Int32(cand), xidx, r, mat, binom)
                return true # break as soon as one exists
            end
        end
    end
    return false
end

# writes out XOR of x and seeds into out 
function xor_many!(st::RState, out::Vector{RankT}, Rn::SetDict,
                   xslot::SlotT, seeds::Vector{SlotT})
    empty!(out)
    k = length(seeds)

    buf = st.buf
    start = 1
    # different cases of when to use xor or xor2
    if xslot != 0
        if k >= 1
            sset_xor2!(out, sd_set(Rn, xslot), sd_set(Rn, @inbounds seeds[1]))
            start = 2
        else
            sx = sd_set(Rn, xslot)
            resize!(out, length(sx)); copyto!(out, sx)
            return out
        end
    elseif k >= 2
        sset_xor2!(out, sd_set(Rn, @inbounds seeds[1]), sd_set(Rn, @inbounds seeds[2]))
        start = 3
    end
    @inbounds for t in start:k
        sset_xor!(out, sd_set(Rn, seeds[t]), buf)
    end
    return out
end

# process boundary of top simplex = xslot U seeds
# returns true if death simplex, false if birth simplex
function process_seeded!(st::RState, xslot::SlotT, seeds::Vector{SlotT}, top_diam::ValueT)::Bool
    n = st.n
    top_bar = st.topbar
    xor_many!(st, top_bar, st.R[n], xslot, seeds)
    isempty(top_bar) && return false
    resolve_bar!(st, n + 1, top_bar, top_diam)
    return true
end

# process if step 2 finds an extra simplex
function process_leftover!(st::RState, xslot::SlotT, top_diam::ValueT)
    n = st.n
    top_bar = st.topbar
    empty!(top_bar)
    append!(top_bar, sd_set(st.R[n], xslot))
    resolve_bar!(st, n + 1, top_bar, top_diam)
    return nothing
end

# run active enumeration on a simplex with xrank, xidx, verts, and sdiam
function run_active_enumeration!(st::RState, xrank::RankT, xidx::CombT,
                                 verts::Vector{Int32}, sdiam::ValueT)
    n = st.n
    m = n + 1
    advance_rbits!(st, sdiam)
    collect_candidates!(st, verts) # step 1 extra vertex candinates
    s = st.seeds
    mat = st.dist.mat
    binom = st.binom
    Rn = st.R[n]
    xslot = sd_slot(Rn, xrank)
    dirty = false
    @inbounds for t in eachindex(s.verts) # for each candidate vertex
        w = s.verts[t]
        sdw = s.seedlists[Int(w) + 1]
        seeded_is_noop(Rn, xslot, sdw) && continue
        if coface_ready!(st.coface_buf, verts, m, w, xidx, sdiam, mat, binom)
            dirty |= process_seeded!(st, xslot, sdw, sdiam)
        end
    end

    # fast path for if step 1 did nothing
    if !dirty
        if has_leftover_coface!(st, verts, m, xidx, sdiam)
            sd_delete!(Rn, xrank)
            sd_delete!(st.Ri[n], xrank)
            return nothing
        end
        register_active!(st.act, xrank, xslot, verts, binom)
        return nothing
    end

    xslot = sd_slot(Rn, xrank) # x could have became inactive via step 1
    if xslot != 0
        if has_leftover_coface!(st, verts, m, xidx, sdiam) # step 2 if still active
            process_leftover!(st, xslot, sdiam)
        end
        xslot = sd_slot(Rn, xrank)
    end

    if xslot != 0
        register_active!(st.act, xrank, xslot, verts, binom) # register active if still active
    end
    return nothing
end

# ============================================================
# Main driver
# ============================================================

function redzed(dist::DistanceMatrix, n::Int; threshold::Real = Inf,
                simplex_pairs::Bool = false)
    n >= 0 || throw(ArgumentError("requires n >= 0"))
    nv = dist.n
    mat = dist.mat

    cone_val = ValueT(Inf) # enclosing radius
    @inbounds for j in 1:nv
        mx = ValueT(-Inf)
        for i in 1:nv
            d = mat[i, j]
            d > mx && (mx = d)
        end
        mx < cone_val && (cone_val = mx)
    end
    thr = min(ValueT(threshold), ValueT(cone_val))

    binom = BinomialTable(nv, n + 2)
    sortsc = SortScratch()

    edges = build_edges(dist, thr, sortsc, binom)

    words = cld(nv, 64)

    face_universe = Int(choose(binom, nv, n))   # (n-1)-simplices have n vertices

    # initialize state 
    st = RState(nv, n, dist, binom,                                    
                [SetDict() for _ in 1:n], [SetDict() for _ in 1:n],    # R, Ri
                [Tuple{ValueT,ValueT}[] for _ in 0:n],                 # pairs
                edges,                                                 # edges
                FaceSlot(0), false, 0,                             # csmap..csdim
                UnionFind(nv),                                         # uf
                RankT[], RankT[], RankT[], RankT[], RankT[], RankT[], RankT[],
                                     # bar buf touched others removed topbar swapbuf
                ActiveReg(n + 1, face_universe),                       # act
                SeedState(nv, words),                                  # seeds
                zeros(UInt64, words, nv),                              # rbits
                fill(Int32(words + 1), nv), zeros(Int32, nv), 1, words, # rlo rhi rptr words
                Vector{Int32}(undef, n + 2), Vector{Int32}(undef, n + 1), # coface_buf top_buf
                [Int32[] for _ in 0:n], [Int32[] for _ in 0:n], simplex_pairs) # sbirth sdeath want_spx

    # dimension 1
    want_h1 = n >= 1
    R1 = want_h1 ? st.R[1] : SetDict()
    Ri1 = want_h1 ? st.Ri[1] : SetDict()
    uf = st.uf
    pairs0 = st.pairs[1]
    ev = edges.verts
    ed = edges.diam
    tb = st.top_buf
    wsx = st.want_spx
    sb0 = st.sbirth[1]
    sd0 = st.sdeath[1]
    @inbounds for e in 1:edges.count
        lo = Int(ev[2 * e - 1])
        hi = Int(ev[2 * e])
        a = uf_find(uf, Int32(lo + 1))
        b = uf_find(uf, Int32(hi + 1))
        ediam = ed[e]
        if a == b
            # birth in dimension 1
            if want_h1
                er = RankT(e)
                sd_singleton!(R1, er, er)
                sd_singleton!(Ri1, er, er)

                # Active enumeration if H1
                if n == 1
                    tb[1] = Int32(lo); tb[2] = Int32(hi)
                    xidx = comb_index(ev, 2 * e - 2, 2, binom)
                    run_active_enumeration!(st, er, xidx, tb, ediam)
                end
            end
        else
            # death in dimension 0: the younger representative dies
            j = max(a, b); o = min(a, b)
            uf.parent[j] = o
            if ediam != 0
                push!(pairs0, (zero(ValueT), ediam))
                if wsx # the younger root is the vertex whose class dies here
                    push!(sb0, j - Int32(1))
                    push!(sd0, Int32(lo), Int32(hi))
                end
            end
        end
    end

    # dimension > 1
    if n >= 2
        for d in 2:n

            uni = Int(choose(binom, nv, d))
            cs = FaceSlot(uni)
            livebuf = RankT[]
            sd_keys!(livebuf, st.R[d - 1])
            @inbounds for z in livebuf
                fs_set!(cs, lvl_idx_of(st, d - 1, z), sd_slot(st.R[d - 1], z))
            end
            st.csmap = cs
            st.cson = true
            st.csdim = d - 1

            m = d + 1
            Rd = st.R[d]
            Rid = st.Ri[d]
            Rprev = st.R[d - 1]
            bar = st.bar
            buf = st.buf
            Rd.side_vpc = m
            Rid.side_vpc = m
            nextrank = Ref(RankT(0))

            # streaming generation one edge at a time 
            generate_stream(st.edges, nv, d, binom,
                function (sdiam::ValueT, cv::Vector{Int32}, cbase::Int)
                    nextrank[] += RankT(1)
                    sr = nextrank[]
                    empty!(bar)
                    @inbounds for omit in 1:m
                        f = face_index_without(cv, cbase, m, omit, binom)
                        fs = fs_get(cs, f)
                        fs != 0 && sset_xor!(bar, sd_set(Rprev, fs), buf)
                    end
                    # on a death this leaves the death simplex in top_buf for resolve_bar!
                    if wsx || isempty(bar)
                        @inbounds for k in 1:m
                            tb[k] = cv[cbase + k]
                        end
                    end
                    if isempty(bar)
                        slot1 = sd_singleton!(Rd, sr, sr)
                        slot2 = sd_singleton!(Rid, sr, sr)
                        sd_set_side!(Rd, slot1, sdiam, tb, 0)
                        sd_set_side!(Rid, slot2, sdiam, tb, 0)
                        if d == n
                            xidx = comb_index(cv, cbase, m, binom)
                            run_active_enumeration!(st, sr, xidx, tb, sdiam)
                        end
                    else
                        resolve_bar!(st, d, bar, sdiam)
                    end
                    return nothing
                end)

            # dimension d-1 final
            st.cson = false
            st.csmap = FaceSlot(0)
            emit_infinite!(st, d - 1)
            st.R[d - 1] = SetDict()
            st.Ri[d - 1] = SetDict()
        end
    end

    # dimension 0 infinite intervals
    @inbounds for v in 1:nv
        if uf_find(uf, Int32(v)) == Int32(v)
            push!(st.pairs[1], (zero(ValueT), ValueT(Inf)))
            if wsx
                push!(sb0, Int32(v - 1))
                push!(sd0, Int32(-1), Int32(-1))
            end
        end
    end
    emit_infinite!(st, n) #  dimension n infinite intervals

    intervals = Dict(d => st.pairs[d + 1] for d in 0:n)
    simplex_pairs || return intervals
    return (intervals, Dict(d => spx_pairs(st, d + 1) for d in 0:n))
end

function emit_infinite!(st::RState, p::Int)
    p < 1 && return nothing
    live = RankT[]
    sd_keys!(live, st.Ri[p])
    sort!(live)
    @inbounds for x in live
        push!(st.pairs[p + 1], (lvl_diam_of(st, p, x), ValueT(Inf)))
        if st.want_spx
            append_verts!(st.sbirth[p + 1], st, p, x)
            for _ in 1:(p + 2) # no death simplex, a full width of -1 marks it infinite
                push!(st.sdeath[p + 1], Int32(-1))
            end
        end
    end
    return nothing
end

# materialize the flat vertex buffers of pairs[k] into (birth_spx, death_spx) tuples
# BW = k vertices per birth simplex (dimension k-1), DW = k+1 per death simplex (dimension k)
function _spx_pairs(sb::Vector{Int32}, sd::Vector{Int32}, ::Val{BW}, ::Val{DW}) where {BW,DW}
    cnt = length(sb) ÷ BW
    out = Vector{Tuple{NTuple{BW,Int32},Union{NTuple{DW,Int32},Tuple{}}}}(undef, cnt)
    @inbounds for t in 1:cnt
        bo = (t - 1) * BW
        do_ = (t - 1) * DW
        b = ntuple(i -> sb[bo + i], Val(BW))
        # infinite intervals have no death simplex, marked by the -1 sentinel
        d = sd[do_ + 1] == Int32(-1) ? () : ntuple(i -> sd[do_ + i], Val(DW))
        out[t] = (b, d)
    end
    return out
end

spx_pairs(st::RState, k::Int) = _spx_pairs(st.sbirth[k], st.sdeath[k], Val(k), Val(k + 1))

function redzed(dist::AbstractMatrix, n::Int; threshold::Real = Inf,
                simplex_pairs::Bool = false)
    return redzed(DistanceMatrix(dist), n; threshold = threshold,
                  simplex_pairs = simplex_pairs)
end



