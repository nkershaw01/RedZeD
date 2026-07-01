#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using ValueT = double;
using IndexT = std::int64_t;

static constexpr ValueT INF_VALUE = std::numeric_limits<ValueT>::infinity();


struct MatrixView {
    const std::vector<std::vector<ValueT>>* nested = nullptr;
    const std::vector<ValueT>* flat = nullptr;
    int n = 0;

    explicit MatrixView(const std::vector<std::vector<ValueT>>& dist)
        : nested(&dist),
          n(validate_nested(dist)) {
    }

    MatrixView(const std::vector<ValueT>& dist, int nverts)
        : flat(&dist),
          n(validate_flat(dist, nverts)) {
    }

    inline ValueT operator()(int i, int j) const {
        if (nested != nullptr) {
            return (*nested)[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }

        return (*flat)[static_cast<std::size_t>(i) * static_cast<std::size_t>(n)
                     + static_cast<std::size_t>(j)];
    }

private:
    static int validate_nested(const std::vector<std::vector<ValueT>>& dist) {
        const int nverts = static_cast<int>(dist.size());

        if (nverts <= 0) {
            throw std::invalid_argument("distance matrix must be nonempty");
        }

        for (int i = 0; i < nverts; ++i) {
            if (static_cast<int>(dist[static_cast<std::size_t>(i)].size()) != nverts) {
                throw std::invalid_argument("distance matrix must be square");
            }
        }

        return nverts;
    }

    static int validate_flat(const std::vector<ValueT>& dist, int nverts) {
        if (nverts <= 0) {
            throw std::invalid_argument("nverts must be positive");
        }

        const std::size_t expected =
            static_cast<std::size_t>(nverts) * static_cast<std::size_t>(nverts);

        if (dist.size() != expected) {
            throw std::invalid_argument(
                "flat distance matrix must have exactly nverts*nverts entries");
        }

        return nverts;
    }
};

// ============================================================
// Binomial table / combinadics
// ============================================================

struct BinomialTable {
    std::vector<IndexT> data;
    int nmax = 0;
    int kmax = 0;

    BinomialTable() = default;

    BinomialTable(int nmax_, int kmax_)
        : data(static_cast<std::size_t>(kmax_ + 1) * static_cast<std::size_t>(nmax_ + 1), 0),
          nmax(nmax_),
          kmax(kmax_) {
        for (int n = 0; n <= nmax; ++n) {
            at(n, 0) = 1;

            const int upto = std::min(n, kmax);
            for (int k = 1; k <= upto; ++k) {
                if (k == n) {
                    at(n, k) = 1;
                } else {
                    at(n, k) = get(n - 1, k - 1) + get(n - 1, k);
                }
            }
        }
    }

    inline IndexT& at(int n, int k) {
        return data[static_cast<std::size_t>(k) * (nmax + 1) + n];
    }

    inline IndexT get(int n, int k) const {
        if (k < 0 || n < 0 || k > n || k > kmax || n > nmax) {
            return 0;
        }
        return data[static_cast<std::size_t>(k) * (nmax + 1) + n];
    }
};

inline IndexT choose(const BinomialTable& binom, int n, int k) {
    return binom.get(n, k);
}

int get_max_vertex(const BinomialTable& binom, IndexT idx, int k, int top_n) {
    int low = k - 1;
    int high = top_n;

    while (low < high) {
        const int mid = (low + high + 1) >> 1;

        if (choose(binom, mid, k) <= idx) {
            low = mid;
        } else {
            high = mid - 1;
        }
    }

    return low;
}

void decode_simplex(
    std::vector<int>& out,
    IndexT idx,
    int dim,
    int n,
    const BinomialTable& binom
) {
    IndexT x = idx;
    int top_n = n - 1;

    for (int k = dim + 1; k >= 2; --k) {
        const int v = get_max_vertex(binom, x, k, top_n);
        out[static_cast<std::size_t>(k - 1)] = v;
        x -= choose(binom, v, k);
        top_n = v - 1;
    }

    out[0] = static_cast<int>(x);
}

// ============================================================
// Boundary enumerator
// ============================================================

struct BoundaryEnumerator {
    IndexT below = 0;
    IndexT above = 0;
    int j = 0;
    int k = 0;
    IndexT idx = 0;
    int dim = 0;

    void set_simplex(IndexT idx_, int dim_, int n) {
        below = idx_;
        above = 0;
        j = n - 1;
        k = dim_;
        idx = idx_;
        dim = dim_;
    }

    inline bool has_next_face() const {
        return k >= 0;
    }

    IndexT next_face_index(const BinomialTable& binom) {
        const int jj = get_max_vertex(binom, below, k + 1, j);
        const IndexT c1 = choose(binom, jj, k + 1);
        const IndexT face_index = above - c1 + below;

        j = jj;
        below -= c1;
        above += choose(binom, jj, k);
        --k;

        return face_index;
    }
};

// ============================================================
// Filtration enumeration
// ============================================================

struct Simplex {
    int dim = 0;
    IndexT idx = 0;
    ValueT diam = 0.0;
};

inline bool simplex_less(const Simplex& a, const Simplex& b) {
    if (a.diam < b.diam) return true;
    if (a.diam > b.diam) return false;

    if (a.dim < b.dim) return true;
    if (a.dim > b.dim) return false;

    return a.idx < b.idx;
}

inline ValueT child_diameter(
    const std::vector<int>& verts,
    int parent_dim,
    ValueT parent_diam,
    int w,
    const MatrixView& dist
) {
    ValueT diam = parent_diam;
    const int m = parent_dim + 1;

    for (int i = 0; i < m; ++i) {
        const ValueT d = dist(w, verts[static_cast<std::size_t>(i)]);
        if (d > diam) {
            diam = d;
        }
    }

    return diam;
}

struct SiblingStream {
    std::vector<Simplex> items;
    std::size_t pos = 0;

    inline bool empty() const {
        return pos >= items.size();
    }

    inline const Simplex& head() const {
        return items[pos];
    }

    inline void advance() {
        ++pos;
    }
};

struct HeadEntry {
    Simplex simplex;
    int id = 0;
};

struct HeadEntryGreater {
    bool operator()(const HeadEntry& a, const HeadEntry& b) const {
        return simplex_less(b.simplex, a.simplex);
    }
};

struct FiltrationEnumerator {
    const MatrixView& dist;
    int n;
    int max_dim;
    ValueT threshold;
    BinomialTable binom;
    std::vector<int> verts;
    std::vector<SiblingStream> streams;
    std::priority_queue<HeadEntry, std::vector<HeadEntry>, HeadEntryGreater> heap;

    FiltrationEnumerator(
        const MatrixView& dist_,
        int max_dim_,
        ValueT threshold_
    )
        : dist(dist_),
          n(dist_.n),
          max_dim(max_dim_),
          threshold(threshold_),
          binom(dist_.n, max_dim_ + 2),
          verts(static_cast<std::size_t>(max_dim_ + 1), 0) {
        std::vector<Simplex> roots(static_cast<std::size_t>(n));

        for (int v = 0; v < n; ++v) {
            roots[static_cast<std::size_t>(v)] = Simplex{0, static_cast<IndexT>(v), 0.0};
        }

        streams.push_back(SiblingStream{std::move(roots), 0});

        if (!streams[0].empty()) {
            heap.push(HeadEntry{streams[0].head(), 0});
        }
    }

    std::optional<SiblingStream> child_stream(const Simplex& s) {
        if (s.dim >= max_dim) {
            return std::nullopt;
        }

        decode_simplex(verts, s.idx, s.dim, n, binom);

        const int m = s.dim + 1;
        const int lastv = verts[static_cast<std::size_t>(m - 1)];

        std::vector<Simplex> children;
        children.reserve(static_cast<std::size_t>(std::max(0, n - lastv - 1)));

        for (int w = lastv + 1; w <= n - 1; ++w) {
            const ValueT diam = child_diameter(verts, s.dim, s.diam, w, dist);

            if (diam <= threshold) {
                const IndexT child_idx = s.idx + choose(binom, w, m + 1);
                children.push_back(Simplex{s.dim + 1, child_idx, diam});
            }
        }

        if (children.empty()) {
            return std::nullopt;
        }

        std::sort(children.begin(), children.end(), simplex_less);
        return SiblingStream{std::move(children), 0};
    }

    void add_stream(SiblingStream&& st) {
        streams.push_back(std::move(st));
        const int sid = static_cast<int>(streams.size()) - 1;

        if (!streams[static_cast<std::size_t>(sid)].empty()) {
            heap.push(HeadEntry{streams[static_cast<std::size_t>(sid)].head(), sid});
        }
    }

    void advance_head(int sid) {
        auto& st = streams[static_cast<std::size_t>(sid)];
        st.advance();

        if (!st.empty()) {
            heap.push(HeadEntry{st.head(), sid});
        }
    }

    std::optional<Simplex> next_simplex() {
        if (heap.empty()) {
            return std::nullopt;
        }

        const HeadEntry head = heap.top();
        heap.pop();

        const Simplex s = head.simplex;

        advance_head(head.id);

        auto child = child_stream(s);
        if (child.has_value()) {
            add_stream(std::move(*child));
        }

        return s;
    }
};

// ============================================================
// Active enumeration helpers
// ============================================================

ValueT diameter(
    std::vector<int>& buf,
    IndexT idx,
    int dim,
    const MatrixView& dist,
    int n,
    const BinomialTable& binom
) {
    decode_simplex(buf, idx, dim, n, binom);

    const int m = dim + 1;
    ValueT diam = 0.0;

    for (int a = 0; a < m; ++a) {
        const int va = buf[static_cast<std::size_t>(a)];

        for (int b = a + 1; b < m; ++b) {
            const ValueT d = dist(va, buf[static_cast<std::size_t>(b)]);
            if (d > diam) {
                diam = d;
            }
        }
    }

    return diam;
}

inline IndexT face_index_without(
    const std::vector<int>& verts,
    int m,
    int omit,
    const BinomialTable& binom
) {
    IndexT idx = 0;

    for (int pos = 0; pos < omit; ++pos) {
        idx += choose(binom, verts[static_cast<std::size_t>(pos)], pos + 1);
    }

    for (int pos = omit + 1; pos < m; ++pos) {
        idx += choose(binom, verts[static_cast<std::size_t>(pos)], pos);
    }

    return idx;
}

inline bool contains_vertex(
    const std::vector<int>& verts,
    int m,
    int w
) {
    for (int i = 0; i < m; ++i) {
        if (verts[static_cast<std::size_t>(i)] == w) {
            return true;
        }
    }
    return false;
}

void remove_face_entry(
    std::vector<std::pair<IndexT, int>>& lst,
    IndexT idx
) {
    for (std::size_t k = 0; k < lst.size(); ++k) {
        if (lst[k].first == idx) {
            lst[k] = lst.back();
            lst.pop_back();
            return;
        }
    }
}

inline int merge_vertex(
    std::vector<int>& out,
    const std::vector<int>& xverts,
    int m,
    int w
) {
    int i = 0;
    int j = 0;

    while (i < m && xverts[static_cast<std::size_t>(i)] < w) {
        out[static_cast<std::size_t>(j)] = xverts[static_cast<std::size_t>(i)];
        ++i;
        ++j;
    }

    out[static_cast<std::size_t>(j)] = w;
    const int wpos = j;
    ++j;

    while (i < m) {
        out[static_cast<std::size_t>(j)] = xverts[static_cast<std::size_t>(i)];
        ++i;
        ++j;
    }

    return wpos;
}

inline ValueT face_diameter(
    const std::vector<int>& verts,
    int m,
    int omit,
    const MatrixView& dist
) {
    ValueT diam = 0.0;

    for (int a = 0; a < m; ++a) {
        if (a == omit) continue;

        const int va = verts[static_cast<std::size_t>(a)];

        for (int b = a + 1; b < m; ++b) {
            if (b == omit) continue;

            const ValueT d = dist(va, verts[static_cast<std::size_t>(b)]);
            if (d > diam) {
                diam = d;
            }
        }
    }

    return diam;
}

bool coface_ready(
    std::vector<int>& candverts,
    const std::vector<int>& xverts,
    int w,
    IndexT xidx,
    ValueT xdiam,
    const MatrixView& dist,
    const BinomialTable& binom
) {
    const int m = static_cast<int>(xverts.size());
    const int total = m + 1;

    const int wpos = merge_vertex(candverts, xverts, m, w);

    for (int omit = 0; omit < total; ++omit) {
        if (omit == wpos) {
            continue;
        }

        const IndexT faceidx = face_index_without(candverts, total, omit, binom);

        if (faceidx < xidx) {
            continue;
        }

        const ValueT facedia = face_diameter(candverts, total, omit, dist);
        if (!(facedia < xdiam)) {
            return false;
        }
    }

    return true;
}

// ============================================================
// Dictionaries
// ============================================================

using MaxSet = std::unordered_set<IndexT>;
using Dict = std::unordered_map<IndexT, MaxSet>;
using Interval = std::pair<ValueT, ValueT>;
using RedzedResult = std::vector<std::vector<Interval>>;

inline void toggle(MaxSet& s, IndexT x) {
    auto it = s.find(x);
    if (it == s.end()) {
        s.insert(x);
    } else {
        s.erase(it);
    }
}

inline void xor_add(MaxSet& dst, const MaxSet& src) {
    for (const IndexT x : src) {
        toggle(dst, x);
    }
}

inline bool xor_from(
    MaxSet& dst,
    const Dict& R,
    IndexT x
) {
    auto it = R.find(x);
    if (it == R.end()) {
        return false;
    }

    xor_add(dst, it->second);
    return true;
}

inline void insert_live(
    Dict& R,
    Dict& Ri,
    IndexT x
) {
    MaxSet sx;
    sx.reserve(1);
    sx.insert(x);
    R[x] = std::move(sx);

    MaxSet rx;
    rx.reserve(1);
    rx.insert(x);
    Ri[x] = std::move(rx);
}

inline bool is_live(
    const Dict& R,
    IndexT x
) {
    return R.find(x) != R.end();
}

void remove_maximal(
    Dict& R,
    Dict& Ri,
    IndexT y,
    std::vector<IndexT>& removed
) {
    removed.clear();

    auto users_it = Ri.find(y);

    if (users_it == Ri.end()) {
        auto set_y_it = R.find(y);

        if (set_y_it != R.end() &&
            set_y_it->second.size() == 1 &&
            set_y_it->second.find(y) != set_y_it->second.end()) {
            R.erase(set_y_it);
            removed.push_back(y);
        }

        return;
    }

    std::vector<IndexT> touched;
    touched.reserve(users_it->second.size());

    for (const IndexT z : users_it->second) {
        touched.push_back(z);
    }

    Ri.erase(users_it);

    for (const IndexT z : touched) {
        auto set_z_it = R.find(z);
        if (set_z_it == R.end()) {
            continue;
        }

        set_z_it->second.erase(y);

        if (set_z_it->second.empty()) {
            R.erase(set_z_it);
            removed.push_back(z);
        }
    }
}

void replace_pivot(
    Dict& R,
    Dict& Ri,
    const MaxSet& bar,
    IndexT j,
    std::vector<IndexT>& removed
) {
    removed.clear();

    auto users_it = Ri.find(j);
    if (users_it == Ri.end()) {
        return;
    }

    std::vector<IndexT> touched;
    touched.reserve(users_it->second.size());

    for (const IndexT z : users_it->second) {
        touched.push_back(z);
    }

    std::vector<IndexT> others;
    others.reserve(bar.size() > 0 ? bar.size() - 1 : 0);

    for (const IndexT y : bar) {
        if (y != j) {
            others.push_back(y);
        }
    }

    for (const IndexT z : touched) {
        auto set_z_it = R.find(z);
        if (set_z_it == R.end()) {
            continue;
        }

        xor_add(set_z_it->second, bar);

        for (const IndexT y : others) {
            auto& set_y = Ri[y];
            toggle(set_y, z);

            if (set_y.empty()) {
                Ri.erase(y);
            }
        }

        if (set_z_it->second.empty()) {
            R.erase(set_z_it);
            removed.push_back(z);
        }
    }

    Ri.erase(j);
}

std::pair<IndexT, ValueT> pick_pivot(
    const MaxSet& bar,
    int dim,
    const MatrixView& dist,
    int nverts,
    const BinomialTable& binom,
    std::vector<int>& buf
) {
    bool first_seen = true;
    IndexT best_idx = 0;
    ValueT best_diam = 0.0;

    for (const IndexT j : bar) {
        const ValueT dj = diameter(buf, j, dim, dist, nverts, binom);

        if (first_seen ||
            dj > best_diam ||
            (dj == best_diam && j > best_idx)) {
            best_idx = j;
            best_diam = dj;
            first_seen = false;
        }
    }

    return {best_idx, best_diam};
}

// ============================================================
// Active enumeration maps
// ============================================================

void register_active(
    std::unordered_map<IndexT, std::vector<int>>& active_verts,
    std::unordered_map<IndexT, std::vector<std::pair<IndexT, int>>>& face_to_active,
    IndexT idx,
    const std::vector<int>& verts,
    int n,
    const BinomialTable& binom
) {
    std::vector<int> stored = verts;
    active_verts[idx] = stored;

    const int m = n + 1;

    for (int omit = 0; omit < m; ++omit) {
        const IndexT faceidx = face_index_without(stored, m, omit, binom);
        face_to_active[faceidx].push_back({idx, stored[static_cast<std::size_t>(omit)]});
    }
}

void unregister_active(
    std::unordered_map<IndexT, std::vector<int>>& active_verts,
    std::unordered_map<IndexT, std::vector<std::pair<IndexT, int>>>& face_to_active,
    const std::vector<IndexT>& doomed,
    int n,
    const BinomialTable& binom
) {
    const int m = n + 1;

    for (const IndexT idx : doomed) {
        auto verts_it = active_verts.find(idx);
        if (verts_it == active_verts.end()) {
            continue;
        }

        const auto& verts = verts_it->second;

        for (int omit = 0; omit < m; ++omit) {
            const IndexT faceidx = face_index_without(verts, m, omit, binom);
            auto lst_it = face_to_active.find(faceidx);

            if (lst_it != face_to_active.end()) {
                remove_face_entry(lst_it->second, idx);

                if (lst_it->second.empty()) {
                    face_to_active.erase(lst_it);
                }
            }
        }

        active_verts.erase(verts_it);
    }
}

// ============================================================
// Seed state
// ============================================================

struct SeedState {
    std::vector<unsigned char> seen;
    std::vector<int> verts;
    std::vector<std::vector<IndexT>> seeds;

    explicit SeedState(int nverts)
        : seen(static_cast<std::size_t>(nverts), 0),
          seeds(static_cast<std::size_t>(nverts)) {
    }

    void reset() {
        for (const int w : verts) {
            seen[static_cast<std::size_t>(w)] = 0;
            seeds[static_cast<std::size_t>(w)].clear();
        }
        verts.clear();
    }
};

void collect_candidates(
    SeedState& state,
    const std::vector<int>& xverts,
    ValueT r,
    const MatrixView& dist,
    int n,
    const BinomialTable& binom,
    const std::unordered_map<IndexT, std::vector<std::pair<IndexT, int>>>& face_to_active
) {
    state.reset();

    const int m = n + 1;

    for (int omit = 0; omit < m; ++omit) {
        const IndexT faceidx = face_index_without(xverts, m, omit, binom);
        auto lst_it = face_to_active.find(faceidx);

        if (lst_it == face_to_active.end()) {
            continue;
        }

        const int xextra = xverts[static_cast<std::size_t>(omit)];

        for (const auto& pair : lst_it->second) {
            const IndexT yidx = pair.first;
            const int w = pair.second;

            if (dist(xextra, w) <= r) {
                if (!state.seen[static_cast<std::size_t>(w)]) {
                    state.seen[static_cast<std::size_t>(w)] = 1;
                    state.verts.push_back(w);
                }

                state.seeds[static_cast<std::size_t>(w)].push_back(yidx);
            }
        }
    }
}

// ============================================================
// Lazy r-close neighbor cache
// ============================================================

struct NeighborCache {
    std::vector<std::vector<int>> order;
    std::vector<int> limit;
};

NeighborCache build_neighbor_cache(
    const MatrixView& dist,
    int nverts
) {
    NeighborCache cache;
    cache.order.resize(static_cast<std::size_t>(nverts));
    cache.limit.assign(static_cast<std::size_t>(nverts), 0);

    for (int v = 0; v < nverts; ++v) {
        auto& row = cache.order[static_cast<std::size_t>(v)];
        row.reserve(static_cast<std::size_t>(std::max(0, nverts - 1)));

        for (int w = 0; w < nverts; ++w) {
            if (w != v) {
                row.push_back(w);
            }
        }

        std::stable_sort(
            row.begin(),
            row.end(),
            [&](int a, int b) {
                return dist(v, a) < dist(v, b);
            });
    }

    return cache;
}

inline int advance_limit(
    NeighborCache& cache,
    int v,
    ValueT r,
    const MatrixView& dist
) {
    auto& row = cache.order[static_cast<std::size_t>(v)];
    int limit = cache.limit[static_cast<std::size_t>(v)];
    const int L = static_cast<int>(row.size());

    while (limit < L && dist(v, row[static_cast<std::size_t>(limit)]) <= r) {
        ++limit;
    }

    cache.limit[static_cast<std::size_t>(v)] = limit;
    return limit;
}

inline std::pair<int, int> choose_pivot_vertex(
    const std::vector<int>& xverts,
    int m,
    ValueT r,
    const MatrixView& dist,
    NeighborCache& cache
) {
    int best_pos = 0;
    int best_lim = std::numeric_limits<int>::max();

    for (int i = 0; i < m; ++i) {
        const int limit = advance_limit(
            cache,
            xverts[static_cast<std::size_t>(i)],
            r,
            dist);

        if (limit < best_lim) {
            best_lim = limit;
            best_pos = i;
        }
    }

    return {best_pos, best_lim};
}

bool has_leftover_coface(
    const std::vector<unsigned char>& seen_w,
    const std::vector<int>& xverts,
    IndexT xidx,
    ValueT r,
    const MatrixView& dist,
    const BinomialTable& binom,
    NeighborCache& neighbors,
    std::vector<int>& candverts
) {
    const int m = static_cast<int>(xverts.size());

    const auto [pivot_pos, limit] =
        choose_pivot_vertex(xverts, m, r, dist, neighbors);

    const int pivot = xverts[static_cast<std::size_t>(pivot_pos)];
    const auto& row = neighbors.order[static_cast<std::size_t>(pivot)];

    for (int t = 0; t < limit; ++t) {
        const int w = row[static_cast<std::size_t>(t)];

        if (contains_vertex(xverts, m, w)) {
            continue;
        }

        if (seen_w[static_cast<std::size_t>(w)]) {
            continue;
        }

        bool ok = true;

        for (int i = 0; i < m; ++i) {
            if (i == pivot_pos) {
                continue;
            }

            if (dist(w, xverts[static_cast<std::size_t>(i)]) > r) {
                ok = false;
                break;
            }
        }

        if (!ok) {
            continue;
        }

        if (coface_ready(candverts, xverts, w, xidx, r, dist, binom)) {
            return true;
        }
    }

    return false;
}

// ============================================================
// Top-level coface processing
// ============================================================

void process_seeded(
    Dict& R,
    Dict& Ri,
    std::vector<Interval>& top_pairs,
    IndexT xidx,
    const std::vector<IndexT>& other_seed_faces,
    ValueT top_diam,
    const MatrixView& dist,
    int nverts,
    const BinomialTable& binom,
    MaxSet& top_bar,
    std::vector<int>& diam_buf,
    std::vector<IndexT>& removed,
    std::unordered_map<IndexT, std::vector<int>>& active_verts,
    std::unordered_map<IndexT, std::vector<std::pair<IndexT, int>>>& face_to_active,
    int n
) {
    top_bar.clear();

    xor_from(top_bar, R, xidx);

    for (const IndexT yidx : other_seed_faces) {
        xor_from(top_bar, R, yidx);
    }

    if (top_bar.empty()) {
        return;
    }

    if (top_bar.size() == 1) {
        const IndexT j = *top_bar.begin();

        remove_maximal(R, Ri, j, removed);
        unregister_active(active_verts, face_to_active, removed, n, binom);

        const ValueT jdiam = diameter(diam_buf, j, n, dist, nverts, binom);

        if (top_diam != jdiam) {
            top_pairs.push_back({jdiam, top_diam});
        }

        return;
    }

    const auto [j, jdiam] = pick_pivot(top_bar, n, dist, nverts, binom, diam_buf);

    replace_pivot(R, Ri, top_bar, j, removed);

    if (!removed.empty()) {
        unregister_active(active_verts, face_to_active, removed, n, binom);
    }

    if (top_diam != jdiam) {
        top_pairs.push_back({jdiam, top_diam});
    }
}

void process_leftover(
    Dict& R,
    Dict& Ri,
    std::vector<Interval>& top_pairs,
    IndexT xidx,
    ValueT top_diam,
    const MatrixView& dist,
    int nverts,
    const BinomialTable& binom,
    std::vector<int>& diam_buf,
    std::vector<IndexT>& removed,
    std::unordered_map<IndexT, std::vector<int>>& active_verts,
    std::unordered_map<IndexT, std::vector<std::pair<IndexT, int>>>& face_to_active,
    int n
) {
    auto r_it = R.find(xidx);
    if (r_it == R.end()) {
        return;
    }

    if (Ri.find(xidx) != Ri.end()) {
        const IndexT j = xidx;

        remove_maximal(R, Ri, j, removed);
        unregister_active(active_verts, face_to_active, removed, n, binom);

        const ValueT jdiam = diameter(diam_buf, j, n, dist, nverts, binom);

        if (top_diam != jdiam) {
            top_pairs.push_back({jdiam, top_diam});
        }

        return;
    }

    MaxSet top_bar = r_it->second;

    if (top_bar.empty()) {
        return;
    }

    if (top_bar.size() == 1) {
        const IndexT j = *top_bar.begin();

        remove_maximal(R, Ri, j, removed);
        unregister_active(active_verts, face_to_active, removed, n, binom);

        const ValueT jdiam = diameter(diam_buf, j, n, dist, nverts, binom);

        if (top_diam != jdiam) {
            top_pairs.push_back({jdiam, top_diam});
        }

        return;
    }

    const auto [j, jdiam] = pick_pivot(top_bar, n, dist, nverts, binom, diam_buf);

    replace_pivot(R, Ri, top_bar, j, removed);

    if (!removed.empty()) {
        unregister_active(active_verts, face_to_active, removed, n, binom);
    }

    if (top_diam != jdiam) {
        top_pairs.push_back({jdiam, top_diam});
    }
}

// ============================================================
// Infinite intervals
// ============================================================

void append_infinite(
    std::vector<Interval>& out,
    const Dict& Ri,
    int dim,
    const MatrixView& dist,
    int nverts,
    const BinomialTable& binom,
    std::vector<int>& diam_buf
) {
    std::vector<IndexT> live_maxima;
    live_maxima.reserve(Ri.size());

    for (const auto& kv : Ri) {
        live_maxima.push_back(kv.first);
    }

    std::sort(live_maxima.begin(), live_maxima.end());

    for (const IndexT x : live_maxima) {
        const ValueT birth = diameter(diam_buf, x, dim, dist, nverts, binom);
        out.push_back({birth, INF_VALUE});
    }
}

// ============================================================
// Main algorithm
// ============================================================

RedzedResult redzed_impl(
    const MatrixView& dist,
    int max_dim,
    ValueT threshold
) {
    if (max_dim < 1) {
        throw std::invalid_argument("redzed requires max_dim >= 1");
    }

    // cone_val = minimum(maximum(dist, dims=1))
    ValueT cone_val = INF_VALUE;

    for (int col = 0; col < dist.n; ++col) {
        ValueT col_max = 0.0;

        for (int row = 0; row < dist.n; ++row) {
            const ValueT d = dist(row, col);
            if (d > col_max) {
                col_max = d;
            }
        }

        if (col_max < cone_val) {
            cone_val = col_max;
        }
    }

    if (cone_val < threshold) {
        threshold = cone_val;
    }

    FiltrationEnumerator fe(dist, max_dim, threshold);
    BoundaryEnumerator be;

    std::vector<Dict> R(static_cast<std::size_t>(max_dim + 1));
    std::vector<Dict> Ri(static_cast<std::size_t>(max_dim + 1));
    RedzedResult pairs(static_cast<std::size_t>(max_dim + 1));

    MaxSet bar;
    MaxSet top_bar;

    std::vector<std::vector<int>> diam_bufs;
    diam_bufs.reserve(static_cast<std::size_t>(max_dim + 1));

    for (int d = 0; d <= max_dim; ++d) {
        diam_bufs.emplace_back(static_cast<std::size_t>(d + 1), 0);
    }

    std::vector<int> top_buf(static_cast<std::size_t>(max_dim + 1), 0);
    std::vector<int> coface_buf(static_cast<std::size_t>(max_dim + 2), 0);

    std::unordered_map<IndexT, std::vector<int>> active_verts;
    std::unordered_map<IndexT, std::vector<std::pair<IndexT, int>>> face_to_active;

    SeedState seeds(fe.n);
    NeighborCache neighbors = build_neighbor_cache(fe.dist, fe.n);

    std::vector<IndexT> removed;

    while (true) {
        auto maybe_s = fe.next_simplex();
        if (!maybe_s.has_value()) {
            break;
        }

        const Simplex s = *maybe_s;

        if (s.dim == 0) {
            insert_live(R[0], Ri[0], s.idx);
            continue;
        }

        bar.clear();
        be.set_simplex(s.idx, s.dim, fe.n);

        while (be.has_next_face()) {
            const IndexT face = be.next_face_index(fe.binom);
            xor_from(bar, R[static_cast<std::size_t>(s.dim - 1)], face);
        }

        if (bar.empty()) {
            const IndexT xidx = s.idx;

            insert_live(
                R[static_cast<std::size_t>(s.dim)],
                Ri[static_cast<std::size_t>(s.dim)],
                xidx);

            if (s.dim == max_dim) {
                decode_simplex(top_buf, s.idx, max_dim, fe.n, fe.binom);

                collect_candidates(
                    seeds,
                    top_buf,
                    s.diam,
                    fe.dist,
                    max_dim,
                    fe.binom,
                    face_to_active);

                for (const int w : seeds.verts) {
                    if (coface_ready(
                            coface_buf,
                            top_buf,
                            w,
                            xidx,
                            s.diam,
                            fe.dist,
                            fe.binom)) {
                        process_seeded(
                            R[static_cast<std::size_t>(max_dim)],
                            Ri[static_cast<std::size_t>(max_dim)],
                            pairs[static_cast<std::size_t>(max_dim)],
                            xidx,
                            seeds.seeds[static_cast<std::size_t>(w)],
                            s.diam,
                            fe.dist,
                            fe.n,
                            fe.binom,
                            top_bar,
                            diam_bufs[static_cast<std::size_t>(max_dim)],
                            removed,
                            active_verts,
                            face_to_active,
                            max_dim);
                    }
                }

                if (is_live(R[static_cast<std::size_t>(max_dim)], xidx)) {
                    if (has_leftover_coface(
                            seeds.seen,
                            top_buf,
                            xidx,
                            s.diam,
                            fe.dist,
                            fe.binom,
                            neighbors,
                            coface_buf)) {
                        process_leftover(
                            R[static_cast<std::size_t>(max_dim)],
                            Ri[static_cast<std::size_t>(max_dim)],
                            pairs[static_cast<std::size_t>(max_dim)],
                            xidx,
                            s.diam,
                            fe.dist,
                            fe.n,
                            fe.binom,
                            diam_bufs[static_cast<std::size_t>(max_dim)],
                            removed,
                            active_verts,
                            face_to_active,
                            max_dim);
                    }
                }

                if (is_live(R[static_cast<std::size_t>(max_dim)], xidx)) {
                    register_active(
                        active_verts,
                        face_to_active,
                        xidx,
                        top_buf,
                        max_dim,
                        fe.binom);
                }
            }

        } else if (bar.size() == 1) {
            const IndexT j = *bar.begin();

            remove_maximal(
                R[static_cast<std::size_t>(s.dim - 1)],
                Ri[static_cast<std::size_t>(s.dim - 1)],
                j,
                removed);

            const ValueT jdiam = diameter(
                diam_bufs[static_cast<std::size_t>(s.dim - 1)],
                j,
                s.dim - 1,
                fe.dist,
                fe.n,
                fe.binom);

            if (s.diam != jdiam) {
                pairs[static_cast<std::size_t>(s.dim - 1)].push_back({jdiam, s.diam});
            }

        } else {
            const auto [j, jdiam] = pick_pivot(
                bar,
                s.dim - 1,
                fe.dist,
                fe.n,
                fe.binom,
                diam_bufs[static_cast<std::size_t>(s.dim - 1)]);

            replace_pivot(
                R[static_cast<std::size_t>(s.dim - 1)],
                Ri[static_cast<std::size_t>(s.dim - 1)],
                bar,
                j,
                removed);

            if (s.diam != jdiam) {
                pairs[static_cast<std::size_t>(s.dim - 1)].push_back({jdiam, s.diam});
            }
        }
    }

    for (int p = 0; p <= max_dim; ++p) {
        append_infinite(
            pairs[static_cast<std::size_t>(p)],
            Ri[static_cast<std::size_t>(p)],
            p,
            fe.dist,
            fe.n,
            fe.binom,
            diam_bufs[static_cast<std::size_t>(p)]);
    }

    return pairs;
}


// ============================================================
// Public API
// ============================================================

RedzedResult redzed(
    const std::vector<std::vector<ValueT>>& dist,
    int max_dim,
    ValueT threshold = INF_VALUE
) {
    MatrixView view(dist);
    return redzed_impl(view, max_dim, threshold);
}

RedzedResult redzed_flat(
    const std::vector<ValueT>& dist,
    int nverts,
    int max_dim,
    ValueT threshold = INF_VALUE
) {
    MatrixView view(dist, nverts);
    return redzed_impl(view, max_dim, threshold);
}
