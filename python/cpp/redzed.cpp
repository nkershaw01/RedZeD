
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <limits>
#include <chrono>
#include <stdexcept>

#ifdef _WIN32
#ifndef NOMINMAX          // libstdc++ already sets it, msvc's stl does not
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
#define REDZED_FSEEK _fseeki64
#define REDZED_FTELL _ftelli64
#else
#define REDZED_FSEEK fseeko
#define REDZED_FTELL ftello
#endif

// count trailing zeros of a nonzero 64 bit word
// __builtin_ctzll is a gcc and clang extension, msvc spells it _BitScanForward64
#if defined(_MSC_VER)
#include <intrin.h>
static inline int redzed_ctz64(uint64_t x) {
    unsigned long i;
    _BitScanForward64(&i, x);
    return (int)i;
}
#else
static inline int redzed_ctz64(uint64_t x) { return __builtin_ctzll(x); }
#endif

typedef float value_t;
typedef int64_t comb_t;    // combinatorial index of a simplex

typedef int64_t rank_t;    // position of a simplex in its dimension's filtration order

typedef int64_t slot_t;

static const value_t VINF = std::numeric_limits<value_t>::infinity();

// (birth_spx, death_spx) as vertex lists, the death list is empty for an infinite interval
typedef std::pair<std::vector<int32_t>, std::vector<int32_t>> SimplexPair;

// ------------------------------------------------------------
// FastMap: rank_t key (must be >= 1) -> slot_t value
// Key -1 marks an empty slot
// Used by SetDict for R and Ri, by the active registry, and by FaceSlot when sparse
// ------------------------------------------------------------

// hand written rather than unordered_map: chaining means a node allocation per entry
// and a pointer chase per lookup, the scattered access this whole design avoids
// shift comes from the capacity, not hard coded: a fixed >> 40 gives only 24 distinct
// hashes and collapses into linear scanning once the table passes 2^24 slots
struct FastMap {
    std::vector<rank_t> keys;
    std::vector<slot_t> vals;
    size_t mask = 15, count = 0, limit = 14;
    int shift = 60;      // 64 - log2(capacity): Fibonacci hashing takes the top bits
    FastMap(size_t cap = 16) { init(cap); }
    void init(size_t cap) {
        size_t c = 16;
        while (c < cap) c <<= 1;
        keys.assign(c, -1); vals.assign(c, 0);
        mask = c - 1; count = 0; limit = (c * 7) >> 3;
        shift = shift_for(c);
    }
    static inline int shift_for(size_t c) {
        int b = 0;
        while (((size_t)1 << b) < c) ++b;
        return 64 - b;
    }
    static inline size_t hsh(rank_t k, int shift) {
        return (size_t)(((uint64_t)(int64_t)k * 0x9E3779B97F4A7C15ull) >> shift); // Fibonacci hashing
    }
    // Get slot for key k
    inline slot_t get(rank_t k) const {
        size_t i = hsh(k, shift);
        for (;;) {
            rank_t kk = keys[i];
            if (kk == k) return vals[i];
            if (kk == -1) return 0;
            i = (i + 1) & mask;
        }
    }
    // Grow when too full
    void grow() {
        std::vector<rank_t> ok = keys;
        std::vector<slot_t> ov = vals;
        size_t c = (mask + 1) << 1;
        keys.assign(c, -1); vals.assign(c, 0);
        mask = c - 1; limit = (c * 7) >> 3; shift = shift_for(c);
        for (size_t t = 0; t < ok.size(); ++t) {
            if (ok[t] == -1) continue;
            size_t i = hsh(ok[t], shift);
            while (keys[i] != -1) i = (i + 1) & mask;
            keys[i] = ok[t]; vals[i] = ov[t];
        }
    }
    // Set a value
    inline void set(rank_t k, slot_t v) {
        size_t i = hsh(k, shift);
        for (;;) {
            rank_t kk = keys[i];
            if (kk == k) { vals[i] = v; return; }
            if (kk == -1) {
                keys[i] = k; vals[i] = v;
                if (++count > limit) grow();
                return;
            }
            i = (i + 1) & mask;
        }
    }
    // returns k's value if present else 0, plus the index k belongs at
    inline slot_t find_or_hole(rank_t k, size_t& hole) const {
        size_t i = hsh(k, shift);
        for (;;) {
            rank_t kk = keys[i];
            if (kk == k) return vals[i];
            if (kk == -1) { hole = i; return 0; }
            i = (i + 1) & mask;
        }
    }
    // Insert at a hole from find_or_hole. Writes before growing, so the hole stays valid
    inline void insert_at(size_t hole, rank_t k, slot_t v) {
        keys[hole] = k; vals[hole] = v;
        if (++count > limit) grow();
    }
    // Backward-shift deletion: keeps the table tombstone-free so lookups never degrade after a long run of births and deaths
    void erase(rank_t k) {
        size_t i = hsh(k, shift);
        for (;;) {
            rank_t kk = keys[i];
            if (kk == -1) return;
            if (kk == k) break;
            i = (i + 1) & mask;
        }
        size_t j = i;
        for (;;) {
            j = (j + 1) & mask;
            rank_t kj = keys[j];
            if (kj == -1) break;
            size_t h = hsh(kj, shift);
            if (((j - h) & mask) >= ((j - i) & mask)) {
                keys[i] = kj; vals[i] = vals[j]; i = j;
            }
        }
        keys[i] = -1; --count;
    }
};

// ------------------------------------------------------------
// Sorted-vector set operations (SSet, ascending)
// ------------------------------------------------------------
typedef std::vector<rank_t> SSet;

// remove x from v
static inline void sset_remove(SSet& v, rank_t x) {
    auto it = std::lower_bound(v.begin(), v.end(), x);
    if (it != v.end() && *it == x) v.erase(it);
}
// toggle x in v
static inline void sset_toggle(SSet& v, rank_t x) {
    auto it = std::lower_bound(v.begin(), v.end(), x);
    if (it != v.end() && *it == x) v.erase(it); else v.insert(it, x);
}
// sorted set xor but output is a seperate vector so no buffer needs to be used
static void sset_xor2(SSet& out, const SSet& a, const SSet& b) {
    size_t na = a.size(), nb = b.size();
    if (out.size() < na + nb) out.resize(na + nb);
    size_t i = 0, j = 0, o = 0;
    while (i < na && j < nb) {
        rank_t x = a[i], y = b[j];
        if (x < y) out[o++] = x, ++i;
        else if (x > y) out[o++] = y, ++j;
        else ++i, ++j;
    }
    while (i < na) out[o++] = a[i++];
    while (j < nb) out[o++] = b[j++];
    out.resize(o);
}
// dst = dst XOR src using sorted sets
// add the smaller element and advance that set, on a tie advance both and add nothing
static void sset_xor(SSet& dst, const SSet& src, SSet& buf) {
    size_t ni = dst.size(), nj = src.size();
    if (nj == 0) return;
    if (ni == 0) { dst = src; return; }
    if (buf.size() < ni + nj) buf.resize(ni + nj);
    size_t i = 0, j = 0, o = 0;
    while (i < ni && j < nj) {
        rank_t a = dst[i], b = src[j];
        if (a < b) buf[o++] = a, ++i;
        else if (a > b) buf[o++] = b, ++j;
        else ++i, ++j;
    }
    while (i < ni) buf[o++] = dst[i++];
    while (j < nj) buf[o++] = src[j++];
    dst.assign(buf.begin(), buf.begin() + o);
}

// ------------------------------------------------------------
// SetDict: rank_t key -> sorted SSet
// ------------------------------------------------------------

// FastMap: Key --> slot
// sets[slot] is corresponding value
// freelist stores slots not in use, if nonempty, pop a slot to store a new value (if empty grow sets)

// sing is vect, sing[i]=0 means sets[i] is not a singleton. Otherwise sing[i]=x means sets[i]={x}
// Benefit is sing is generally small enough to store in cache
struct SetDict {
    FastMap map;
    std::vector<SSet> sets;
    std::vector<slot_t> freelist;
    std::vector<rank_t> sing;

    int side_vpc = 0;
    std::vector<value_t> sdiam;
    std::vector<int32_t> svert;
    // record a key's own diameter and vertex list against its pool slot
    inline void set_side(slot_t s, value_t d, const int32_t* v) {
        if (!side_vpc) return;
        sdiam[s - 1] = d;
        std::memcpy(&svert[(size_t)(s - 1) * side_vpc], v,
                    sizeof(int32_t) * side_vpc);
    }
    // rank -> its stored diameter / vertex list for live keys only
    inline value_t side_diam(rank_t k) const { return sdiam[map.get(k) - 1]; }
    inline const int32_t* side_vert(rank_t k) const {
        return &svert[(size_t)(map.get(k) - 1) * side_vpc];
    }
    SetDict() : map(64) {}

    inline slot_t slot(rank_t k) const { return map.get(k); }       // key -> slot
    inline SSet& at(slot_t s) { return sets[s - 1]; }               // slot -> sets[slot]
    inline const SSet& at(slot_t s) const { return sets[s - 1]; }
    inline bool haskey(rank_t k) const { return map.get(k) != 0; }  // Check if key is a key
    // update sing after sets[slot] has changed
    inline void touch(slot_t s) {
        const SSet& v = sets[s - 1];
        sing[s - 1] = v.size() == 1 ? v[0] : 0;
    }
    // Creates (or clears) the entry for k and returns its slot
    slot_t create(rank_t k) {
        size_t hole = 0;
        slot_t s = map.find_or_hole(k, hole);
        if (s) { sets[s - 1].clear(); sing[s - 1] = 0; return s; } // defensive, k already a key, clear its slot and return
        if (freelist.empty()) {
            sets.push_back(SSet());
            sing.push_back(0);
            if (side_vpc) { sdiam.push_back(0);
                            svert.resize(sets.size() * (size_t)side_vpc); }
            s = (slot_t)sets.size();
        } else {
            s = freelist.back(); freelist.pop_back();
            sets[s - 1].clear();
        }
        sing[s - 1] = 0;
        map.insert_at(hole, k, s);
        return s;
    }
    // Creates the entry k = {v}
    slot_t singleton(rank_t k, rank_t v) {
        slot_t s = create(k);       // create the vector for k, keeping its slot
        sets[s - 1].push_back(v);   // add v
        sing[s - 1] = v;            // update sing
        return s;
    }
    // delete k from dict
    void erase(rank_t k) {
        slot_t s = map.get(k);
        if (!s) return;
        map.erase(k);
        sets[s - 1].clear();
        sing[s - 1] = 0;
        freelist.push_back(s);
    }
    // remove y from sets[slot]
    inline void remove_elem(slot_t s, rank_t y) {
        sset_remove(sets[s - 1], y);
        touch(s);
    }
    // sets[slot] <- sets[slot] XOR src, by merging into scratch and then
    // swapping the two vectors rather than copying the result back
    inline void xor_into(slot_t s, const SSet& src, SSet& scratch) {
        SSet& dst = sets[s - 1];
        if (src.empty()) return;
        if (dst.empty()) { dst = src; touch(s); return; }
        sset_xor2(scratch, dst, src);
        dst.swap(scratch);
        touch(s);
    }
    // Iterates the live keys
    void keys(std::vector<rank_t>& out) const {
        out.clear();
        for (size_t i = 0; i < map.keys.size(); ++i)
            if (map.keys[i] != -1) out.push_back(map.keys[i]);
    }
};

// fast check done during active enumeration to avoid the entire xor if it would become zero (even if a face is active)
// observation is that the majority of the time this happens, it is something like {a}xor{a} with singletons
// this check determines if you would just get cancellations via singletons
// it is fast since sing is generally small and fast to access
// xslot = slot of x, the simplex active enumeration is run on
// beneficial to run before coface_ready to avoid the expensive check that x is the youngest face when possible too

static inline bool seeded_is_noop(const SetDict& sd, slot_t xslot,
                                  const std::vector<slot_t>& seeds) {
    int cnt = 0;
    rank_t f = 0;
    if (xslot) {
        f = sd.sing[xslot - 1];
        if (!f) return false;   // not a singleton
        cnt = 1;
    }
    for (size_t t = 0; t < seeds.size(); ++t) {
        rank_t v = sd.sing[seeds[t] - 1];
        if (!v) return false;
        if (cnt == 0) f = v;
        else if (v != f) return false;
        ++cnt;
    }
    return (cnt & 1) == 0;
}

// ------------------------------------------------------------
// UnionFind for dimension 0
//
// R_0[v] is always the singleton {rep(v)} and Ri_0[r] is
// r's connected component

// equivalent to union find where root is forced to be older
// ------------------------------------------------------------

struct UnionFind {
    std::vector<int32_t> parent;    // parent[v] is v's parent
    explicit UnionFind(int n) : parent(n) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    // find the root of x
    inline int32_t find(int32_t x) {
        int32_t r = x;
        while (parent[r] != r) r = parent[r];        // first pass finding root
        while (parent[x] != r) { int32_t nx = parent[x]; parent[x] = r; x = nx; } // second pass updating path
        return r;
    }
};

// ------------------------------------------------------------
// FaceSlot: combinatorial index -> reduction-set pool slot (0 = not present)
// faces come as vertex sets, so this gets R[face] by going
// comb index -> slot -> pool[slot] = R[face]

// also used by active enumeration, given a face gets
// {set of all active n-simps contianin that face} via
// comb index -> slot ->
// ar.lr[slot]  = ranks of actives containing it
// ar.lw[slot]  = each one's extra vertex
// ar.ls[slot]  = each one's R-pool slot

// when dense gets too large, switches to sparse
// ------------------------------------------------------------

static const int64_t DENSE_INDEX_CAP = 1 << 24;

// ------------------------------------------------------------
// FaceSlot: combinatorial index -> pool slot (0 = not present)
// ------------------------------------------------------------

struct FaceSlot {
    std::vector<slot_t> dense;
    FastMap sparse;
    bool isdense = true;
    FaceSlot() {}

    void init(int64_t universe) {
        dense.clear(); dense.shrink_to_fit();
        sparse.init(16);
        if (universe >= 0 && universe <= DENSE_INDEX_CAP) {
            isdense = true;
            dense.assign((size_t)universe, 0);
        } else {
            isdense = false;
        }
    }
    inline slot_t get(comb_t k) const {
        if (isdense) return dense[(size_t)k];
        return sparse.get(k);
    }
    inline void set(comb_t k, slot_t v) {
        if (isdense) dense[(size_t)k] = v;
        else if (v == 0) sparse.erase(k);   // 0 = delete
        else sparse.set(k, v);
    }
};

// ------------------------------------------------------------
// Distance matrix
// ------------------------------------------------------------

struct DistanceMatrix {
    int n = 0;      // num verts
    std::vector<value_t> mat;
    inline void resize_lower(int nv) { n = nv; mat.resize((size_t)nv * (nv - 1) / 2); }
    static inline size_t offset(int i, int j) { return (size_t)j * (j - 1) / 2 + i; }
    inline value_t operator()(int i, int j) const {
        if (i == j) return (value_t)0;
        return i < j ? mat[offset(i, j)] : mat[offset(j, i)];
    }
};

// ------------------------------------------------------------
// Binomial table
// ------------------------------------------------------------

struct BinomialTable {
    std::vector<comb_t> data;
    int kmax = 0;
    BinomialTable() {}
    BinomialTable(int nmax, int kmaxv) : kmax(kmaxv) {
        data.assign((size_t)(kmax + 1) * (nmax + 1), 0);
        for (int n = 0; n <= nmax; ++n) {
            data[(size_t)n * (kmax + 1)] = 1;
            for (int k = 1; k <= std::min(n, kmax); ++k) {
                if (k == n) data[(size_t)n * (kmax + 1) + k] = 1;
                else data[(size_t)n * (kmax + 1) + k] =
                         data[(size_t)(n - 1) * (kmax + 1) + (k - 1)] +
                         data[(size_t)(n - 1) * (kmax + 1) + k];
            }
        }
    }
    inline comb_t operator()(int n, int k) const {
        if (k > n || k < 0 || n < 0) return 0;
        return data[(size_t)n * (kmax + 1) + k];
    }
};

// Combinatorial index of an ascending vertex list of length m
static inline comb_t comb_index(const int32_t* v, int m, const BinomialTable& b) {
    comb_t idx = 0;
    for (int k = 1; k <= m; ++k) idx += b(v[k - 1], k);
    return idx;
}
// Combinatorial index of the face obtained by omitting position omit
static inline comb_t face_index_without(const int32_t* v, int m, int omit,
                                        const BinomialTable& b) {
    comb_t idx = 0;
    for (int k = 1; k <= omit - 1; ++k) idx += b(v[k - 1], k);
    for (int k = omit + 1; k <= m; ++k) idx += b(v[k - 1], k - 1);
    return idx;
}

// ------------------------------------------------------------
// EdgeList: the dimension-1 simplices in filtration order
// ------------------------------------------------------------

struct EdgeList {
    size_t count = 0;
    std::vector<value_t> diam;      // diameters
    std::vector<int32_t> verts;     // vertices, flat, 2 per edge: edge r is [2r-2, 2r-1]
    inline const int32_t* vp(size_t r) const { return &verts[2 * (r - 1)]; }
    inline value_t dm(size_t r) const { return diam[r - 1]; }
};

// ------------------------------------------------------------
// Sorting a level into filtration order
// ------------------------------------------------------------

struct SortScratch {
    std::vector<rank_t> cnt;    // histogram of bucket counts
    std::vector<value_t> dbuf, sdbuf;
    std::vector<int32_t> vbuf, svbuf;
    std::vector<std::pair<comb_t, rank_t>> kv;
    std::vector<value_t> kvd;
};

// ------------------------------------------------------------
// Bucket sort of a generated level into filtration order
// ------------------------------------------------------------

// find the bucket based on diameter d. lo = range min, inv = 1/(hi-lo), nb = num buckets
static inline int bucket_of(value_t d, value_t lo, value_t inv, int nb) {
    int b = (int)((d - lo) * inv * nb);
    if (b < 0) return 0;
    if (b >= nb) return nb - 1;
    return b;
}

// direct sort for n<=32, depth >=5, or if dim_dim==dim_max
// comparison sort used when n>32, else insertion sort
static void sort_slice_small(std::vector<value_t>& diam, std::vector<int32_t>& verts, int vpc,
                             size_t a, size_t b, const BinomialTable& bt, SortScratch& sc,
                             bool comparison) {
    size_t n = b - a + 1;
    if (n <= 1) return;
    if (sc.kv.size() < n) { sc.kv.resize(n); sc.kvd.resize(n); }
    if (sc.sdbuf.size() < n) sc.sdbuf.resize(n);
    if (sc.svbuf.size() < n * (size_t)vpc) sc.svbuf.resize(n * (size_t)vpc);
    for (size_t t = 0; t < n; ++t) {
        size_t p = a + t - 1;
        sc.kvd[t] = diam[p];
        sc.kv[t] = { comb_index(&verts[p * vpc], vpc, bt), (rank_t)p };
    }
    std::vector<rank_t>& ord = sc.cnt;
    if (ord.size() < n) ord.resize(n);
    for (size_t t = 0; t < n; ++t) ord[t] = (rank_t)t;
    auto less = [&](rank_t x, rank_t y) {
        if (sc.kvd[x] != sc.kvd[y]) return sc.kvd[x] < sc.kvd[y];
        return sc.kv[x].first < sc.kv[y].first;
    };
    if (comparison) std::sort(ord.begin(), ord.begin() + n, less);  // built in comparison sort
    else {
        // insertion sort
        for (size_t i = 1; i < n; ++i) {
            rank_t x = ord[i];
            size_t j = i;
            while (j > 0 && less(x, ord[j - 1])) { ord[j] = ord[j - 1]; --j; }
            ord[j] = x;
        }
    }
    // order in buffer
    for (size_t t = 0; t < n; ++t) {
        size_t p = (size_t)sc.kv[ord[t]].second;
        sc.sdbuf[t] = diam[p];
        std::memcpy(&sc.svbuf[t * vpc], &verts[p * vpc], sizeof(int32_t) * vpc);
    }
    // copy back
    for (size_t t = 0; t < n; ++t) {
        diam[a + t - 1] = sc.sdbuf[t];
        std::memcpy(&verts[(a + t - 1) * vpc], &sc.svbuf[t * vpc], sizeof(int32_t) * vpc);
    }
}

// sort a bucket, a = start index, b = end index
static void sort_slice(std::vector<value_t>& diam, std::vector<int32_t>& verts, int vpc,
                       size_t a, size_t b, const BinomialTable& bt, SortScratch& sc,
                       int depth) {
    size_t n = b - a + 1;
    if (n <= 1) return;     // singleton
    if (n <= 32 || depth >= 5) { sort_slice_small(diam, verts, vpc, a, b, bt, sc, n > 32); return; } // small or depth >=5 directly sort

    value_t dmin = VINF, dmax = -VINF;
    for (size_t i = a; i <= b; ++i) {
        value_t d = diam[i - 1];
        if (d < dmin) dmin = d;
        if (d > dmax) dmax = d;
    }
    if (dmax == dmin) { sort_slice_small(diam, verts, vpc, a, b, bt, sc, true); return; } // recursive on diam offers nothing

    int nb = (int)(n >> 2);     // n/4 buckets
    if (nb < 2) nb = 2;         // floor of two, unreachable but safe if n cap changes
    value_t inv = 1.0 / (dmax - dmin);
    std::vector<rank_t> cnt(nb + 1, 0);
    for (size_t i = a; i <= b; ++i) cnt[bucket_of(diam[i - 1], dmin, inv, nb) + 1] += 1;   // make histogram
    for (int t = 1; t <= nb; ++t) cnt[t] += cnt[t - 1];  // cnt[0] = 0, cnt[bk] = number of elements in buckets 0 â€¦ bk-1 = bucket bk's 0-based start offset
    std::vector<rank_t> bounds(cnt);

    // resize diam and vertex buffs
    if (sc.dbuf.size() < n) sc.dbuf.resize(n);
    if (sc.vbuf.size() < n * (size_t)vpc) sc.vbuf.resize(n * (size_t)vpc);
    // sort into buckets
    for (size_t i = a; i <= b; ++i) {
        value_t d = diam[i - 1];
        int bk = bucket_of(d, dmin, inv, nb);
        rank_t p = cnt[bk]++;
        sc.dbuf[p] = d;
        std::memcpy(&sc.vbuf[(size_t)p * vpc], &verts[(i - 1) * vpc], sizeof(int32_t) * vpc);
    }
    // copy back
    for (size_t t = 0; t < n; ++t) {
        diam[a + t - 1] = sc.dbuf[t];
        std::memcpy(&verts[(a + t - 1) * vpc], &sc.vbuf[t * vpc], sizeof(int32_t) * vpc);
    }
    // recursion on each new bucket
    for (int bk = 0; bk < nb; ++bk) {
        size_t sa = a + (size_t)bounds[bk];
        size_t sb = a + (size_t)bounds[bk + 1] - 1;
        if (sb > sa) sort_slice(diam, verts, vpc, sa, sb, bt, sc, depth + 1);
    }
}

// diam = diameters
// verts = vertices
// vpc = vertices per simplex
// starts = starts of buckets
// nb = number of buckets
static void finish_buckets(std::vector<value_t>& diam, std::vector<int32_t>& verts, int vpc,
                           const std::vector<rank_t>& starts, int nb,
                           const BinomialTable& bt, SortScratch& sc) {
    for (int b = 0; b < nb; ++b) {
        size_t a = (size_t)starts[b] + 1;
        size_t e = (size_t)starts[b + 1];
        if (e > a) sort_slice(diam, verts, vpc, a, e, bt, sc, 1);
    }
}

// ------------------------------------------------------------
// Dimension 1: edges, in filtration order
// first scan builds the diameter histogram, second puts edges in buckets
// ------------------------------------------------------------

static void build_edges(const DistanceMatrix& dist, value_t threshold,
                        const BinomialTable& bt, SortScratch& sc, EdgeList& out) {
    int n = dist.n;
    int64_t cap = (int64_t)n * (n - 1) / 2;
    int nb = (int)(cap >> 2);           // 4 edges per bucket
    if (nb < 1024) nb = 1024;           // lower bound
    if (nb > (1 << 16)) nb = 1 << 16;   // upper bound for memory
    value_t dmin = 0.0;
    value_t inv = (std::isfinite(threshold) && threshold > 0) ? 1.0 / threshold : 0.0;
    std::vector<rank_t> cnt(nb + 1, 0);

    // pass 1: diameter histogram
    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < j; ++i) {
            value_t d = dist(i, j);
            if (d <= threshold) {
                cnt[bucket_of(d, dmin, inv, nb) + 1] += 1;
            }
        }
    }
    for (int t = 1; t <= nb; ++t) cnt[t] += cnt[t - 1];
    size_t count = (size_t)cnt[nb];
    out.count = count;
    if (count == 0) return;
    std::vector<rank_t> starts(cnt);

    out.diam.resize(count);
    out.verts.resize(2 * count);
    // pass 2
    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < j; ++i) {
            value_t d = dist(i, j);
            if (d <= threshold) {
                int b = bucket_of(d, dmin, inv, nb);
                rank_t p = cnt[b]++;
                out.diam[p] = d;
                out.verts[2 * (size_t)p] = i;
                out.verts[2 * (size_t)p + 1] = j;
            }
        }
    }
    finish_buckets(out.diam, out.verts, 2, starts, nb, bt, sc);
}

// ------------------------------------------------------------
// Dimensions >= 2: edge-driven streaming generation
// ------------------------------------------------------------

// g = matrix storing G_<e, the graph with all edges less than current edge
// need = number of vertices needed to form a d simplex when union with edge
// emit = function to process simplex
template <class F>
struct EdStream {
    const uint64_t* g = nullptr;
    int words = 0, n = 0, vpc = 0;
    const BinomialTable* bt = nullptr;
    value_t diam = 0;
    int32_t base[2] = {0, 0};
    std::vector<int32_t> stack, merged;
    std::vector<uint64_t> scratch;
    bool accumulate = false;
    std::vector<comb_t> akey;
    std::vector<int32_t> avert;
    F* emit = nullptr;

    // merge {u,v} with stack[0..k-1] (descending) into ascending
    inline void put(int k) {
        int32_t* m = merged.data();
        int i = 0, j = k - 1, t = 0;
        while (i < 2 && j >= 0)
            m[t++] = (base[i] < stack[j]) ? base[i++] : stack[j--];
        while (i < 2) m[t++] = base[i++];
        while (j >= 0) m[t++] = stack[j--];
        if (accumulate) {   // multiple edges have same length, store comb index to process them with streams created by other edges
            akey.push_back(comb_index(m, vpc, *bt));
            avert.insert(avert.end(), m, m + vpc);
        } else {            // edge is alone at that length
            (*emit)(diam, m);
        }
    }

    void rec(uint64_t* cand, int need, int depth) {
        for (int w = 0; w < words; ++w) {   // number of words
            // this depth's candidate set, word by word
            // candidate set = \cap of N(x) for all current verts, but only those verts less than the last vert added
            uint64_t wb = cand[w];
            while (wb) {
                int bp = redzed_ctz64(wb);      // lowest bit position
                wb &= wb - 1;                   // clear that position
                int x = w * 64 + bp;            // vertex idx
                if (x >= n) break;              // safeguard against leftover higher bits, nothing ever sets those bits so never actually fired
                stack[depth] = (int32_t)x;      // this depths vert
                if (need == 1) { put(depth + 1); }
                else {
                    uint64_t* nxt = &scratch[(size_t)(depth + 1) * words];
                    const uint64_t* gx = g + (size_t)x * words;
                    for (int t = 0; t < words; ++t) nxt[t] = cand[t] & gx[t];   // intersection with x
                    int cw = x >> 6;            // keep only vertices below x
                    for (int t = cw + 1; t < words; ++t) nxt[t] = 0;
                    int rem = x & 63;
                    if (cw < words) nxt[cw] &= (rem ? (((uint64_t)1 << rem) - 1) : 0);
                    rec(nxt, need - 1, depth + 1);
                }
            }
        }
    }

    // sort and process all generated simps from ties
    void flush() {
        size_t c = akey.size();
        if (!c) return;
        std::vector<uint32_t> ord(c);
        for (size_t i = 0; i < c; ++i) ord[i] = (uint32_t)i;
        std::sort(ord.begin(), ord.end(),
                  [&](uint32_t a, uint32_t b) { return akey[a] < akey[b]; });
        for (size_t i = 0; i < c; ++i)
            (*emit)(diam, &avert[(size_t)ord[i] * vpc]);
        akey.clear(); avert.clear();
    }
};

template <class F>
static void generate_stream(const EdgeList& edges, const DistanceMatrix& dist,
                            int d, const BinomialTable& bt, F&& emit) {
    int n = dist.n;
    int words = (n + 63) / 64;
    if (edges.count == 0) return;
    std::vector<uint64_t> g((size_t)words * n, 0);
    int need = d - 1;

    EdStream<F> C;
    C.g = g.data(); C.words = words; C.n = n; C.bt = &bt; C.vpc = d + 1;
    C.stack.assign(need + 1, 0);
    C.merged.assign(d + 2, 0);
    C.scratch.assign((size_t)(need + 1) * words, 0);
    C.emit = &emit;

    size_t E = edges.count, i = 0;
    while (i < E) {
        value_t r = edges.diam[i];
        size_t jend = i;
        // need to walk in edge length tie groups
        while (jend < E && edges.diam[jend] == r) ++jend;
        C.diam = r;
        C.accumulate = (jend != i + 1);     // false if no ties, true otherwise
        for (size_t k = i; k < jend; ++k) {
            int u = edges.verts[2 * k], v = edges.verts[2 * k + 1];
            const uint64_t* gu = g.data() + (size_t)u * words;
            const uint64_t* gv = g.data() + (size_t)v * words;
            uint64_t* cand = &C.scratch[0];
            uint64_t any = 0;
            for (int t = 0; t < words; ++t) { cand[t] = gu[t] & gv[t]; any |= cand[t]; }   // intersection of words
            C.base[0] = (int32_t)u; C.base[1] = (int32_t)v;
            if (any) C.rec(cand, need, 0);
            // update g after edge has been processed
            g[(size_t)u * words + (v >> 6)] |= (uint64_t)1 << (v & 63);
            g[(size_t)v * words + (u >> 6)] |= (uint64_t)1 << (u & 63);
        }
        if (C.accumulate) C.flush();
        i = jend;
    }
}

//############################################################################
// RedZeD_VR algorithm
//############################################################################

// ============================================================
// Active-simplex registry
//
// The keys of R_n are exactly the active n-simplices
// The registry keeps its vertex list
// for every (n-1)-face the list of active simplices containing it (with the extra vertex)
// ============================================================

struct ActiveReg {
    int m = 0;                              // vertices per top simplex (= n + 1)
    FastMap slot_of;                        // active rank -> vertex-slot
    std::vector<int32_t> verts;             // flat, m per slot (VERTEX ids, not ranks)
    std::vector<slot_t> freelist;
    size_t nslots = 0;
    FaceSlot face_slot;                    // face comb index -> list slot
    std::vector<std::vector<rank_t>> lr;    // per list slot: active ranks
    std::vector<std::vector<int32_t>> lw;   // per list slot: the extra vertex
    std::vector<std::vector<slot_t>> ls;    // per list slot: that active's R-set pool slot
    std::vector<slot_t> lfree;
    ActiveReg() : slot_of(1024) {}
};

// ============================================================
// Whole-run state
// ============================================================

struct State {
    int nv = 0, n = 0;                  // vertex count, homology degree
    DistanceMatrix dist;                // distance matrix
    BinomialTable binom;                // binomial table
    std::vector<SetDict> R, Ri;         // vect of R's and Ri's as set dicts
    std::vector<std::vector<std::pair<value_t, value_t>>> pairs;  // pairs found
    EdgeList edges;                     // sorted edges

    // face --> pool slot lookup
    FaceSlot csmap;
    bool cson = false;
    int csdim = 0;

    UnionFind* uf = nullptr;            // union find for dimension 0

    // bar = r \partial x (called bar because notation was origonally r\partial = \overline{\partial})
    // buf = buffer for xor, touched = snapshot of Ri[y], others = bar \ {j} : entries whose Ri needs updating
    // topbar = r \partial for n+1 dimensional simplices, swapbuf = buffer
    SSet bar, buf, touched, others, topbar, swapbuf;
    std::vector<rank_t> removed;        // keys removed from R
    ActiveReg act;

    // Reusable Step-1 state
    std::vector<char> seen;                         // seen verts as a vect
    std::vector<uint64_t> seenbits;                 // same info as a bitset to make step 2 subtraction fast
    std::vector<int32_t> seedverts;                 // set of candidate vertices
    std::vector<std::vector<slot_t>> seedlists;     // R-set slots of the face-sharing actives

    std::vector<uint64_t> rbits;        // bitset adjacency grown one edge at a time
    std::vector<int32_t> rlo, rhi;      // lowest / highest bit for a vertex
    size_t rptr = 1;                    // pointer into edge list
    int words = 0;                      // num words
    // coface_buf = buffer for n+1 simplex, top_buf = buffer for n simplex
    std::vector<int32_t> coface_buf, top_buf;

    // optional simplex pairs: one entry per stored interval of pairs[k] in the same order
    std::vector<std::vector<int32_t>> sbirth;   // flat birth-simplex vertices, k+1 per entry of pairs[k]
    std::vector<std::vector<int32_t>> sdeath;   // flat death-simplex vertices, k+2 per entry, all -1 = infinite
    bool want_spx = false;                      // record simplex pairs, left empty when false
};

// is b ~ a
static inline bool rbit(const State& st, int a, int b) {
    return (st.rbits[(size_t)a * st.words + (b >> 6)] >> (b & 63)) & 1ull;
}

// ============================================================
// Reduction bookkeeping
// ============================================================

// get the diameter of a p dimensional simplex r
// Ri used for diams since only ever called on alive elements
static inline value_t lvl_diam_of(State& st, int p, rank_t r) {
    if (p == 1) return st.edges.dm((size_t)r);
    return st.Ri[p].side_diam(r);
}
// get the vertices of a p dimensional simplex r
static inline const int32_t* lvl_verts_of(State& st, int p, rank_t r) {
    if (p == 1) return st.edges.vp((size_t)r);
    return st.R[p].side_vert(r);
}

// append the vertices of the p dimensional simplex r to out
static inline void append_verts(std::vector<int32_t>& out, State& st, int p, rank_t r) {
    const int32_t* v = (p == 1) ? st.edges.vp((size_t)r) : st.Ri[p].side_vert(r);
    out.insert(out.end(), v, v + (p + 1));
}

// remove z from csmap (comb index --> slot)
// apply sometimes false since there isn't always a csmap (i.e. H1 or during AE)
static inline void cs_clear(State& st, rank_t z, bool apply) {
    if (!apply) return;
    st.csmap.set(comb_index(lvl_verts_of(st, st.csdim, z), st.csdim + 1, st.binom), 0);
}

// special case for if r \partial = {y}, remove y from each R[x] st y in R[y]
static void remove_maximal(State& st, SetDict& Rd, SetDict& Rid, rank_t y, bool apply_cs) {
    st.removed.clear();
    slot_t us = Rid.slot(y);
    st.touched = Rid.at(us);    // Ri[y]
    Rid.erase(y);               // y dies
    for (rank_t z : st.touched) {
        slot_t zs = Rd.slot(z); // never 0: z in Ri[y] means y in R[z], so z is a key of R
        Rd.remove_elem(zs, y);
        if (Rd.at(zs).empty()) {    // z now inactive, remove
            cs_clear(st, z, apply_cs);
            Rd.erase(z);
            st.removed.push_back(z);
        }
    }
}

// general case, bar = r partial size geq 2, j = youngest element
static void replace_pivot(State& st, SetDict& Rd, SetDict& Rid, const SSet& bar, rank_t j,
                          bool apply_cs) {
    st.removed.clear();
    slot_t us = Rid.slot(j);    // never 0: j = bar.back(), same as remove_maximal
    st.touched = Rid.at(us);    // Ri[j]

    st.others.clear();          // others = bar \ {j}
    for (rank_t y : bar) if (y != j) st.others.push_back(y);

    for (rank_t z : st.touched) {
        slot_t zs = Rd.slot(z); // never 0: z in Ri[j] means j in R[z], so z is a key of R
        Rd.xor_into(zs, bar, st.swapbuf);
        if (Rd.at(zs).empty()) {
            cs_clear(st, z, apply_cs);
            Rd.erase(z);
            st.removed.push_back(z);
        }
    }

    size_t K = st.touched.size();
    for (rank_t y : st.others) {
        // never 0: y in bar means y is an element of some R set, so a key of Ri
        SSet& sy = Rid.at(Rid.slot(y));
        if (K * 128 < sy.size()) {  // toggle is faster when touched is significalty smaller than sy
            for (rank_t z : st.touched) sset_toggle(sy, z);
        } else {
            sset_xor(sy, st.touched, st.buf);
        }
        if (sy.empty()) Rid.erase(y);
    }
    Rid.erase(j);
}

// doomed = set of ranks that just became inactive, remove them form active reg if present
static void unregister_active(State& st, const std::vector<rank_t>& doomed) {
    ActiveReg& ar = st.act;
    for (rank_t rk : doomed) {
        slot_t slot = ar.slot_of.get(rk);
        if (!slot) continue;
        const int32_t* base = &ar.verts[(size_t)(slot - 1) * ar.m];
        for (int omit = 1; omit <= ar.m; ++omit) {
            comb_t f = face_index_without(base, ar.m, omit, st.binom);
            slot_t lsx = ar.face_slot.get(f);
            if (!lsx) continue;
            auto& lr = ar.lr[lsx - 1];
            auto& lw = ar.lw[lsx - 1];
            auto& lsl = ar.ls[lsx - 1];
            for (size_t k = 0; k < lr.size(); ++k) {
                if (lr[k] == rk) {
                    lr[k] = lr.back(); lr.pop_back();
                    lw[k] = lw.back(); lw.pop_back();
                    lsl[k] = lsl.back(); lsl.pop_back();
                    break;
                }
            }
            if (lr.empty()) { ar.face_slot.set(f, 0); ar.lfree.push_back(lsx); }
        }
        ar.slot_of.erase(rk);
        ar.freelist.push_back(slot);
    }
}

// register an active simplex that is active after active enumeration
static void register_active(State& st, rank_t rk, slot_t rslot, const int32_t* verts) {
    ActiveReg& ar = st.act;
    slot_t slot;
    // allocate a new slot if freelist is empty, otherwise pop a slot
    if (ar.freelist.empty()) {
        ++ar.nslots;
        ar.verts.resize(ar.nslots * (size_t)ar.m);
        slot = (slot_t)ar.nslots;
    } else {
        slot = ar.freelist.back(); ar.freelist.pop_back();
    }
    int32_t* base = &ar.verts[(size_t)(slot - 1) * ar.m];
    for (int i = 0; i < ar.m; ++i) base[i] = verts[i];
    ar.slot_of.set(rk, slot);
    for (int omit = 1; omit <= ar.m; ++omit) {
        comb_t f = face_index_without(base, ar.m, omit, st.binom);
        slot_t lsx = ar.face_slot.get(f);
        if (!lsx) {
            if (ar.lfree.empty()) {
                ar.lr.push_back({}); ar.lw.push_back({}); ar.ls.push_back({});
                lsx = (slot_t)ar.lr.size();
            } else {
                lsx = ar.lfree.back(); ar.lfree.pop_back();
            }
            ar.face_slot.set(f, lsx);
        }
        ar.lr[lsx - 1].push_back(rk);
        ar.lw[lsx - 1].push_back(base[omit - 1]);
        ar.ls[lsx - 1].push_back(rslot);
    }
}

// handles r partial x \neq 0 for x of dimension d
// sdiam = simplex diam, used to store pair
static void resolve_bar(State& st, int d, SSet& bar, value_t sdiam) {
    SetDict& Rd = st.R[d - 1];
    SetDict& Rid = st.Ri[d - 1];
    bool apply_cs = st.cson && (d - 1 == st.csdim);
    rank_t j = bar.back();
    value_t jdiam = lvl_diam_of(st, d - 1, j);

    // j is the birth simplex of the interval, recorded before its pool slot is freed
    // death simplex is the current d dimensional simplex
    // left in coface_buf by coface_ready in the top dimension and in top_buf by the generator below it
    if (st.want_spx && sdiam != jdiam) {
        append_verts(st.sbirth[d - 1], st, d - 1, j);
        const int32_t* dv = (d == st.n + 1) ? st.coface_buf.data() : st.top_buf.data();
        st.sdeath[d - 1].insert(st.sdeath[d - 1].end(), dv, dv + (d + 1));
    }

    if (bar.size() == 1) remove_maximal(st, Rd, Rid, j, apply_cs);
    else replace_pivot(st, Rd, Rid, bar, j, apply_cs);

    if (d == st.n + 1 && !st.removed.empty()) unregister_active(st, st.removed);  // top dimension, need to unregister no longer active
    if (sdiam != jdiam) st.pairs[d - 1].push_back({jdiam, sdiam});   // store pairs of nonzero length
}

// ============================================================
// Active enumeration
// ============================================================

// advance rbits to new radius r
static void advance_rbits(State& st, value_t r) {
    const EdgeList& lv = st.edges;
    size_t p = st.rptr, cnt = lv.count;
    while (p <= cnt && lv.diam[p - 1] <= r) {
        int i = lv.verts[2 * (p - 1)];
        int j = lv.verts[2 * (p - 1) + 1];
        int wj = j >> 6, wi = i >> 6;
        st.rbits[(size_t)i * st.words + wj] |= (uint64_t)1 << (j & 63);
        st.rbits[(size_t)j * st.words + wi] |= (uint64_t)1 << (i & 63);
        if (st.rlo[i] > wj) st.rlo[i] = wj;
        if (st.rhi[i] < wj) st.rhi[i] = wj;
        if (st.rlo[j] > wi) st.rlo[j] = wi;
        if (st.rhi[j] < wi) st.rhi[j] = wi;
        ++p;
    }
    st.rptr = p;
}

// add a vertex to a sortecd list of vertices
static inline int merge_vertex(int32_t* out, const int32_t* xv, int m, int32_t w) {
    int i = 0, j = 0;
    while (i < m && xv[i] < w) out[j++] = xv[i++];
    out[j] = w;
    int wpos = j + 1;
    ++j;
    while (i < m) out[j++] = xv[i++];
    return wpos;
}

// diameter of the face given by omitting the vertex omit
static inline value_t face_diameter(const int32_t* v, int m, int omit,
                                    const DistanceMatrix& dist) {
    value_t diam = 0.0;
    for (int a = 1; a <= m; ++a) {
        if (a == omit) continue;
        for (int b = a + 1; b <= m; ++b) {
            if (b == omit) continue;
            value_t d = dist(v[a - 1], v[b - 1]);
            if (d > diam) diam = d;
        }
    }
    return diam;
}

// Is xverts U {w} ready to appear, i.e. is xverts its youngest face?
static bool coface_ready(int32_t* cand, const int32_t* xverts, int m, int32_t w, comb_t xidx,
                         value_t xdiam, const DistanceMatrix& dist, const BinomialTable& bt) {
    int total = m + 1;
    int wpos = merge_vertex(cand, xverts, m, w);
    for (int omit = 1; omit <= total; ++omit) {
        if (omit == wpos) continue;
        comb_t f = face_index_without(cand, total, omit, bt);
        // xverts U {w} always has diam = diam xverts, so if f<xidx, face already born
        if (f < xidx) continue;
        if (!(face_diameter(cand, total, omit, dist) < xdiam)) return false;  // if f > xidx, needs to have smaller diam
    }
    return true;
}

// Step 1: candidate vertices w such that xverts U {w} is a valid coface containing
// another active simplex
static void collect_candidates(State& st, const int32_t* xverts) {
    for (int32_t w : st.seedverts) {    // only need to reset vertices actually found, why seedverts is stored
        st.seen[w] = 0;
        st.seedlists[w].clear();
        st.seenbits[w >> 6] = 0;
    }
    st.seedverts.clear();
    ActiveReg& ar = st.act;
    int m = ar.m;
    for (int omit = 1; omit <= m; ++omit) {     // faces of x
        comb_t f = face_index_without(xverts, m, omit, st.binom);
        slot_t lsx = ar.face_slot.get(f);
        if (!lsx) continue;
        int xextra = xverts[omit - 1];
        const auto& lw = ar.lw[lsx - 1];    // extra verts
        const auto& lsl = ar.ls[lsx - 1];   // pool slots
        for (size_t t = 0; t < lw.size(); ++t) {
            int w = lw[t];
            if (rbit(st, xextra, w)) {  // d(xextra,w)<=r, store
                if (!st.seen[w]) {
                    st.seen[w] = 1;
                    st.seedverts.push_back(w);
                    st.seenbits[w >> 6] |= (uint64_t)1 << (w & 63);
                }
                st.seedlists[w].push_back(lsl[t]);  // add (x \ omit U {w}): by end of fcn stores all other active faces of x U w
            }
        }
    }
}

// Step 2: does another coface exist not yet processed
static bool has_leftover_coface(State& st, const int32_t* xverts, int m, comb_t xidx,
                                value_t r) {
    const uint64_t* rb = st.rbits.data();
    int words = st.words, nv = st.nv;
    // only words nonzero for all verts can survive intersection
    int wlo = 0, whi = words - 1;
    for (int i = 0; i < m; ++i) {
        int v = xverts[i];
        if (st.rlo[v] > wlo) wlo = st.rlo[v];
        if (st.rhi[v] < whi) whi = st.rhi[v];
    }
    for (int w = wlo; w <= whi; ++w) {
        uint64_t acc = rb[(size_t)xverts[0] * words + w];
        if (!acc) continue;
        for (int i = 1; i < m; ++i) {   // repeated intersection
            acc &= rb[(size_t)xverts[i] * words + w];
            if (!acc) break;
        }
        if (!acc) continue;
        acc &= ~st.seenbits[w];     // intersect with compliment of seen
        while (acc) {
            int bp = redzed_ctz64(acc);
            acc &= acc - 1;
            int cand = w * 64 + bp;
            if (cand >= nv) break;
            if (coface_ready(st.coface_buf.data(), xverts, m, (int32_t)cand, xidx, r,
                             st.dist, st.binom))
                return true;    // break as soon as one exists
        }
    }
    return false;
}

// writes out XOR of x and seeds into out
static void xor_many(State& st, SSet& out, SetDict& Rn, slot_t xslot,
                     const std::vector<slot_t>& seeds) {
    out.clear();
    size_t k = seeds.size();
    size_t start = 0;
    // different cases of when to use xor or xor2
    if (xslot) {
        if (k >= 1) { sset_xor2(out, Rn.at(xslot), Rn.at(seeds[0])); start = 1; }
        else { out = Rn.at(xslot); return; }
    } else if (k >= 2) {
        sset_xor2(out, Rn.at(seeds[0]), Rn.at(seeds[1])); start = 2;
    }
    for (size_t t = start; t < k; ++t) sset_xor(out, Rn.at(seeds[t]), st.buf);
}

// process boundary of top simplex = xslot U seeds
// returns true if death simplex, false if birth simplex
static bool process_seeded(State& st, slot_t xslot, const std::vector<slot_t>& seeds,
                           value_t top_diam) {
    xor_many(st, st.topbar, st.R[st.n], xslot, seeds);
    if (st.topbar.empty()) return false;
    resolve_bar(st, st.n + 1, st.topbar, top_diam);
    return true;
}

// process if step 2 finds an extra simplex
static void process_leftover(State& st, slot_t xslot, value_t top_diam) {
    st.topbar = st.R[st.n].at(xslot);
    resolve_bar(st, st.n + 1, st.topbar, top_diam);
}

// run active enumeration on a simplex with xrank, xidx, verts, and sdiam
static void run_active_enumeration(State& st, rank_t xrank, comb_t xidx, const int32_t* verts,
                                   value_t sdiam) {
    int n = st.n, m = n + 1;
    advance_rbits(st, sdiam);
    collect_candidates(st, verts);      // step 1 extra vertex candinates
    SetDict& Rn = st.R[n];
    slot_t xslot = Rn.slot(xrank);

    bool dirty = false;
    for (size_t t = 0; t < st.seedverts.size(); ++t) {  // for each candidate vertex
        int32_t w = st.seedverts[t];
        const std::vector<slot_t>& sdw = st.seedlists[w];
        if (seeded_is_noop(Rn, xslot, sdw)) continue;
        if (coface_ready(st.coface_buf.data(), verts, m, w, xidx, sdiam, st.dist, st.binom))
            dirty |= process_seeded(st, xslot, sdw, sdiam);
    }

    // fast path for if step 1 did nothing
    if (!dirty) {
        if (has_leftover_coface(st, verts, m, xidx, sdiam)) {
            Rn.erase(xrank);
            st.Ri[n].erase(xrank);
            return;
        }
        register_active(st, xrank, xslot, verts);
        return;
    }

    xslot = Rn.slot(xrank);     // x could have became inactive via step 1
    if (xslot) {
        if (has_leftover_coface(st, verts, m, xidx, sdiam)) process_leftover(st, xslot, sdiam); // step 2 if still active
        xslot = Rn.slot(xrank);
    }
    if (xslot) register_active(st, xrank, xslot, verts);    // register active if still active
}

static void emit_infinite(State& st, int p) {
    if (p < 1) return;
    std::vector<rank_t> live;
    st.Ri[p].keys(live);
    std::sort(live.begin(), live.end());
    for (rank_t x : live) {
        st.pairs[p].push_back({ lvl_diam_of(st, p, x), VINF });
        if (st.want_spx) {
            append_verts(st.sbirth[p], st, p, x);
            // no death simplex, a full width of -1 marks it infinite
            st.sdeath[p].insert(st.sdeath[p].end(), (size_t)(p + 2), -1);
        }
    }
}

// materialize the flat vertex buffers of pairs[k] into (birth_spx, death_spx) vertex lists
// bw = k+1 vertices per birth simplex (dimension k), dw = k+2 per death simplex (dimension k+1)
static std::vector<SimplexPair> spx_pairs(State& st, int k) {
    const std::vector<int32_t>& sb = st.sbirth[k];
    const std::vector<int32_t>& sd = st.sdeath[k];
    size_t bw = (size_t)k + 1, dw = (size_t)k + 2;
    size_t cnt = sb.size() / bw;
    std::vector<SimplexPair> out(cnt);
    for (size_t t = 0; t < cnt; ++t) {
        out[t].first.assign(sb.begin() + t * bw, sb.begin() + (t + 1) * bw);
        // infinite intervals have no death simplex, marked by the -1 sentinel
        if (sd[t * dw] != -1)
            out[t].second.assign(sd.begin() + t * dw, sd.begin() + (t + 1) * dw);
    }
    return out;
}

// ============================================================
// Main driver
// ============================================================

// out_spx, when not null, receives the simplex pairs, one entry per interval of pairs[d]
static std::vector<std::vector<std::pair<value_t, value_t>>>
redzed(DistanceMatrix& distIn, int n, value_t threshold, bool simplex_pairs = false,
       std::vector<std::vector<SimplexPair>>* out_spx = nullptr) {
    int nv = distIn.n;
    value_t cone = VINF;    // enclosing radius
    for (int j = 0; j < nv; ++j) {
        value_t mx = -VINF;
        for (int i = 0; i < nv; ++i) mx = std::max(mx, distIn(i, j));
        cone = std::min(cone, mx);
    }
    value_t thr = std::min(threshold, cone);

    // initialize state
    State st;
    st.nv = nv; st.n = n;
    st.dist = std::move(distIn);
    st.binom = BinomialTable(nv, n + 2);
    SortScratch sortsc;

    build_edges(st.dist, thr, st.binom, sortsc, st.edges);

    st.words = (nv + 63) / 64;
    st.R.resize(n + 1);
    st.Ri.resize(n + 1);
    st.pairs.assign(n + 1, {});
    st.sbirth.assign(n + 1, {});
    st.sdeath.assign(n + 1, {});
    st.want_spx = simplex_pairs;
    UnionFind uf(nv);
    st.uf = &uf;
    st.act.m = n + 1;
    {
        comb_t fu = st.binom(nv, n);    // (n-1)-simplices have n vertices
        st.act.face_slot.init(fu);
    }
    st.seen.assign(nv, 0);
    st.seenbits.assign(st.words, 0);
    st.seedlists.assign(nv, {});
    st.rbits.assign((size_t)st.words * nv, 0);
    st.rlo.assign(nv, st.words);
    st.rhi.assign(nv, -1);
    st.coface_buf.assign(n + 2, 0);
    st.top_buf.assign(n + 1, 0);

    // dimension 1
    {
        const bool want_h1 = n >= 1;
        SetDict* R1 = want_h1 ? &st.R[1] : nullptr;
        SetDict* Ri1 = want_h1 ? &st.Ri[1] : nullptr;
        const EdgeList& ed = st.edges;
        bool wsx = st.want_spx;
        std::vector<int32_t>& sb0 = st.sbirth[0];
        std::vector<int32_t>& sd0 = st.sdeath[0];
        for (size_t e = 1; e <= ed.count; ++e) {
            int lo = ed.verts[2 * (e - 1)];
            int hi = ed.verts[2 * (e - 1) + 1];
            int32_t a = uf.find(lo), b = uf.find(hi);
            value_t ediam = ed.diam[e - 1];
            if (a == b) {
                // birth in dimension 1
                if (want_h1) {
                    rank_t er = (rank_t)e;
                    R1->singleton(er, er);
                    Ri1->singleton(er, er);

                    // Active enumeration if H1
                    if (n == 1) {
                        st.top_buf[0] = lo; st.top_buf[1] = hi;
                        comb_t xidx = comb_index(&ed.verts[2 * (e - 1)], 2, st.binom);
                        run_active_enumeration(st, er, xidx, st.top_buf.data(), ediam);
                    }
                }
            } else {
                // death in dimension 0: the younger representative dies
                int32_t jj = std::max(a, b), oo = std::min(a, b);
                uf.parent[jj] = oo;
                if (ediam != 0.0) {
                    st.pairs[0].push_back({ 0.0, ediam });
                    if (wsx) {  // the younger root is the vertex whose class dies here
                        sb0.push_back(jj);
                        sd0.push_back((int32_t)lo); sd0.push_back((int32_t)hi);
                    }
                }
            }
        }
    }

    // dimension > 1
    if (n >= 2) {
        for (int d = 2; d <= n; ++d) {
            comb_t uni = st.binom(nv, d);
            st.csmap.init(uni);
            {
                std::vector<rank_t> live;
                st.R[d - 1].keys(live);
                for (rank_t z : live)
                    st.csmap.set(comb_index(lvl_verts_of(st, d - 1, z), d, st.binom),
                                 st.R[d - 1].slot(z));
            }
            st.cson = true;
            st.csdim = d - 1;

            int m = d + 1;
            SetDict& Rd = st.R[d];
            SetDict& Rid = st.Ri[d];
            SetDict& Rprev = st.R[d - 1];
            Rd.side_vpc = m;
            Rid.side_vpc = m;
            rank_t next_rank = 0;

            // streaming generation one edge at a time
            generate_stream(st.edges, st.dist, d, st.binom,
                [&](value_t sdiam, const int32_t* cv) {
                    rank_t sr = ++next_rank;
                    st.bar.clear();
                    for (int omit = 1; omit <= m; ++omit) {
                        comb_t f = face_index_without(cv, m, omit, st.binom);
                        slot_t fs = st.csmap.get(f);
                        if (fs) sset_xor(st.bar, Rprev.at(fs), st.buf);
                    }
                    if (st.bar.empty()) {
                        slot_t s1 = Rd.singleton(sr, sr);
                        slot_t s2 = Rid.singleton(sr, sr);
                        Rd.set_side(s1, sdiam, cv);
                        Rid.set_side(s2, sdiam, cv);
                        if (d == n) {
                            comb_t xidx = comb_index(cv, m, st.binom);
                            run_active_enumeration(st, sr, xidx, cv, sdiam);
                        }
                    } else {
                        // on a death this leaves the death simplex in top_buf for resolve_bar
                        if (st.want_spx)
                            std::memcpy(st.top_buf.data(), cv, sizeof(int32_t) * (size_t)m);
                        resolve_bar(st, d, st.bar, sdiam);
                    }
                });

            // dimension d-1 final
            st.cson = false;
            st.csmap.init(0);
            emit_infinite(st, d - 1);
            st.R[d - 1] = SetDict();
            st.Ri[d - 1] = SetDict();
        }
    }

    // dimension 0 infinite intervals
    for (int v = 0; v < nv; ++v)
        if (uf.find(v) == v) {
            st.pairs[0].push_back({ 0.0, VINF });
            if (st.want_spx) {
                st.sbirth[0].push_back((int32_t)v);
                st.sdeath[0].push_back(-1); st.sdeath[0].push_back(-1);
            }
        }
    emit_infinite(st, n);   // dimension n infinite intervals

    if (simplex_pairs && out_spx) {
        out_spx->assign(n + 1, {});
        for (int d = 0; d <= n; ++d) (*out_spx)[d] = spx_pairs(st, d);
    }
    return st.pairs;
}

// the python extension compiles this file with -DREDZED_NO_MAIN to drop the command
// line driver below: a shared library must not define main. cli builds define nothing
// and are unaffected
#ifndef REDZED_NO_MAIN

// usage: redzed <matrix.bin> <dim> [--lower] [--timing] [--threshold R] [--barcode]
//        [--simplex-pairs]
// threshold defaults to +inf, which redzed then caps at the enclosing radius
int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr,
                "usage: %s <matrix.bin> <dim> [--lower] [--timing] [--threshold R] "
                "[--barcode] [--simplex-pairs]\n",
                argv[0]);
        return 1;
    }
    bool timing = false;
    bool lower_input = false;
    bool simplex_pairs = false;
    bool barcode = false;
    value_t threshold = VINF;
    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--timing") timing = true;
        else if (a == "--simplex-pairs") simplex_pairs = true;
        else if (a == "--barcode") barcode = true;
        else if (a == "--lower") lower_input = true;
        else if (a == "--threshold" && i + 1 < argc) threshold = (value_t)atof(argv[++i]);
        else { fprintf(stderr, "unrecognized argument: %s\n", a.c_str()); return 1; }
    }

    FILE* f = fopen(argv[1], "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    DistanceMatrix D;
    if (lower_input) {
        if (REDZED_FSEEK(f, 0, SEEK_END) != 0) { fprintf(stderr, "seek failed\n"); return 1; }
        int64_t bytes = (int64_t)REDZED_FTELL(f);
        REDZED_FSEEK(f, 0, SEEK_SET);
        int64_t m = bytes / (int64_t)sizeof(value_t);
        int64_t nv = (int64_t)((1.0 + std::sqrt(1.0 + 8.0 * (double)m)) / 2.0 + 0.5);
        if (nv < 1 || nv * (nv - 1) / 2 != m) {
            fprintf(stderr, "%s holds %lld values, which is not a lower triangle\n",
                    argv[1], (long long)m);
            return 1;
        }
        D.resize_lower((int)nv);
        if (m && fread(D.mat.data(), sizeof(value_t), (size_t)m, f) != (size_t)m) return 1;
    } else {
        int64_t n64 = 0;
        if (fread(&n64, sizeof(int64_t), 1, f) != 1) return 1;
        D.resize_lower((int)n64);
        std::vector<double> row((size_t)n64);
        for (int64_t j = 0; j < n64; ++j) {
            if (fread(row.data(), sizeof(double), row.size(), f) != row.size()) return 1;
            value_t* dst = D.mat.data() + (size_t)(j * (j - 1) / 2);
            for (int64_t i = 0; i < j; ++i) dst[i] = (value_t)row[(size_t)i];
        }
    }
    fclose(f);

    int dim = atoi(argv[2]);
    std::vector<std::vector<SimplexPair>> spx;
    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::pair<value_t, value_t>>> pairs =
        redzed(D, dim, threshold, simplex_pairs, &spx);
    auto t1 = std::chrono::high_resolution_clock::now();

    for (int d = 0; d <= dim; ++d) {
        double chk = 0.0;
        long long ninf = 0;
        std::vector<std::pair<value_t, value_t>> v = pairs[d];
        std::sort(v.begin(), v.end());
        for (auto& pr : v) {
            if (std::isinf(pr.second)) { ++ninf; chk += pr.first; }
            else chk += pr.first * 3.0 + pr.second * 7.0;
        }
        printf("H%d %zu %lld %.10f\n", d, v.size(), ninf, chk);
    }
    // BAR <dim> <birth> <death>, in the same order as the SPX lines below
    if (barcode) {
        for (int d = 0; d <= dim; ++d)
            for (size_t t = 0; t < pairs[d].size(); ++t)
                printf("BAR %d %.10g %.10g\n", d, (double)pairs[d][t].first,
                       (double)pairs[d][t].second);
    }
    // SPX <dim> <birth> <death> <birth simplex> | <death simplex>, empty death = infinite
    if (simplex_pairs) {
        for (int d = 0; d <= dim; ++d)
            for (size_t t = 0; t < spx[d].size(); ++t) {
                printf("SPX %d %.10g %.10g", d, (double)pairs[d][t].first,
                       (double)pairs[d][t].second);
                for (int32_t vv : spx[d][t].first) printf(" %d", vv);
                printf(" |");
                for (int32_t vv : spx[d][t].second) printf(" %d", vv);
                printf("\n");
            }
    }
    if (timing) {
        double sec = std::chrono::duration<double>(t1 - t0).count();
        fprintf(stderr, "TIME %.6f\n", sec);
#ifdef _WIN32
        PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
            fprintf(stderr, "PEAKRSS %llu\n", (unsigned long long)pmc.PeakWorkingSetSize);
#endif
    }
    return 0;
}

#endif  // REDZED_NO_MAIN


