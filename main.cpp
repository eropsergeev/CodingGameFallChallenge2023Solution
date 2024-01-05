// #define DETERMINISTIC_MODE

#include <bits/stdc++.h>

using namespace std;

namespace {

template <class T>
struct Vec2 {
    T x, y;
    template <class U, std::enable_if_t<!std::is_same_v<U, std::common_type_t<T, U>>, bool> = true>
    [[gnu::always_inline]] explicit operator Vec2<U>() const noexcept {
        return Vec2<U>{static_cast<U>(x), static_cast<U>(y)};
    }
    template <class U, std::enable_if_t<std::is_same_v<U, std::common_type_t<T, U>>, bool> = true>
    [[gnu::always_inline]] operator Vec2<U>() const noexcept {
        return Vec2<U>{static_cast<U>(x), static_cast<U>(y)};
    }
    [[gnu::always_inline]] constexpr Vec2() = default;
    [[gnu::always_inline]] constexpr Vec2(T x, T y) : x(x), y(y) {}
    [[gnu::always_inline]] explicit Vec2(double radians): x(std::cos(radians)), y(std::sin(radians)) {}
    template <class V>
    [[gnu::always_inline]] explicit constexpr Vec2(const V& v) noexcept : x(v.x), y(v.y) {}
    [[gnu::always_inline]] Vec2& operator+=(const Vec2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }
    [[gnu::always_inline]] Vec2& operator-=(const Vec2& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }
    template <class U>
    [[gnu::always_inline, nodiscard, gnu::pure]] Vec2<std::common_type_t<T, U>> operator+(
        const Vec2<U>& other) const noexcept {
        return {x + other.x, y + other.y};
    }
    template <class U>
    [[gnu::always_inline, nodiscard, gnu::pure]] Vec2<std::common_type_t<T, U>> operator-(
        const Vec2<U>& other) const noexcept {
        return {x - other.x, y - other.y};
    }
    template<class U = std::conditional_t<std::is_floating_point_v<T>, T, double>>
    [[gnu::always_inline, nodiscard, gnu::pure]] Vec2<U> normalize() const noexcept {
        Vec2<U> ret(*this);
        ret /= ret.len();
        return ret;
    }
    [[gnu::always_inline, nodiscard, gnu::pure]] T lenSq() const noexcept { return x * x + y * y; }
    [[gnu::always_inline, nodiscard, gnu::pure]] auto len() const noexcept {
        return std::sqrt(lenSq());
    }
    [[gnu::always_inline, nodiscard, gnu::pure]] Vec2 operator-() const noexcept {
        return Vec2{-x, -y};
    }
    [[gnu::always_inline]] Vec2& operator*=(T c) {
        x *= c;
        y *= c;
        return *this;
    }
    template <class Dummy = T, std::enable_if_t<std::is_floating_point_v<Dummy>, bool> = true>
    [[gnu::always_inline]] Vec2& operator/=(T c) {
        x /= c;
        y /= c;
        return *this;
    }
    template <class U>
    [[gnu::always_inline, nodiscard, gnu::pure]] auto rotate(U radians) const {
        using Ret = std::conditional_t<std::is_floating_point_v<T>, T, double>;
        using Ang = std::common_type_t<Ret, U>;
        auto sn = std::sin(static_cast<Ang>(radians));
        auto cs = std::cos(static_cast<Ang>(radians));
        return Vec2<Ret>{x * cs - y * sn, x * sn + y * cs};
    }
    [[gnu::always_inline, nodiscard, gnu::pure]] Vec2 left() const noexcept {
        Vec2 ans = *this;
        std::swap(ans.x, ans.y);
        ans.x = -ans.x;
        return ans;
    }
    [[gnu::always_inline, nodiscard, gnu::pure]] Vec2 right() const noexcept {
        Vec2 ans = *this;
        std::swap(ans.x, ans.y);
        ans.y = -ans.y;
        return ans;
    }
    [[gnu::always_inline, nodiscard, gnu::pure]] auto radians() const { return std::atan2(y, x); }
    [[gnu::always_inline, nodiscard, gnu::pure]] bool operator==(const Vec2& other) const noexcept {
        return x == other.x && y == other.y;
    }
    [[gnu::always_inline, nodiscard, gnu::pure]] bool operator!=(const Vec2& other) const noexcept {
        return !(*this == other);
    }
};

template <class T1, class T2>
[[gnu::always_inline, nodiscard, gnu::pure]] inline Vec2<std::common_type_t<T1, T2>> operator*(
    const Vec2<T1>& v, T2 c) noexcept {
    return {v.x * c, v.y * c};
}

template <class T1, class T2>
[[gnu::always_inline, nodiscard, gnu::pure]] inline auto operator*(T2 c,
                                                                   const Vec2<T1>& v) noexcept {
    return v * c;
}

template <class T1, class T2>
[[gnu::always_inline, nodiscard, gnu::pure]] inline auto operator/(const Vec2<T1>& v,
                                                                   T2 c) noexcept {
    using Ret = std::conditional_t<std::is_integral_v<T1> && std::is_integral_v<T2>,
                                   double,
                                   std::common_type_t<T1, T2>>;
    return Vec2<Ret>{static_cast<Ret>(v.x) / c, static_cast<Ret>(v.y) / c};
}

template <class T1, class T2>
[[gnu::always_inline, nodiscard, gnu::pure]] inline bool near(const Vec2<T1>& a, const Vec2<T2>& b,
                                                              double eps = 0) noexcept {
    return near(a.x, b.x, eps) && near(a.y, b.y, eps);
}

template <class T1, class T2>
[[gnu::always_inline, nodiscard, gnu::pure]] inline auto dot(const Vec2<T1>& a,
                                                             const Vec2<T2>& b) noexcept {
    return a.x * b.x + a.y * b.y;
}

template <class T1, class T2>
[[gnu::always_inline, nodiscard, gnu::pure]] inline auto operator*(const Vec2<T1>& a,
                                                             const Vec2<T2>& b) noexcept {
    return dot(a, b);
}

template <class T1, class T2>
[[gnu::always_inline, nodiscard, gnu::pure]] inline auto cross(const Vec2<T1>& a,
                                                               const Vec2<T2>& b) noexcept {
    return a.x * b.y - a.y * b.x;
}

template <class T>
inline std::ostream& operator<<(std::ostream& out, const Vec2<T>& v) {
    return out << v.x << " " << v.y;
}

template <class T>
inline std::istream& operator>>(std::istream& in, Vec2<T>& v) {
    return in >> v.x >> v.y;
}

using Vec2f = Vec2<float>;
using Vec2d = Vec2<double>;
using Vec2i = Vec2<int>;
using Vec2i32 = Vec2<int32_t>;
using Vec2i64 = Vec2<int64_t>;

#define REGISTER_CMD(_mnemonic, ...) struct Cmd##_mnemonic : public tuple<__VA_ARGS__> { \
    using Base = tuple<__VA_ARGS__>; \
    static constexpr string_view mnemonic = #_mnemonic; \
}

struct Rect {
    Vec2i p1, p2;
    Rect intersect(const Rect &other) const {
        Rect ret{{max(p1.x, other.p1.x), max(p1.y, other.p1.y)}, {min(p2.x, other.p2.x), min(p2.y, other.p2.y)}};
        if (ret.p1.x > ret.p2.x || ret.p1.y > ret.p2.y)
            return Rect{};
        return ret;
    }
    bool contains(Vec2i p) const {
        return p1.x <= p.x && p.x <= p2.x && p1.y <= p.y && p.y <= p1.y;
    }
    int area() const {
        return width() * height();
    }
    int width() const {
        return p2.x - p1.x;
    }
    int height() const {
        return p2.y - p1.y;
    }
    array<Vec2i, 4> points() const {
        return {p1, {p1.x, p2.y}, p2, {p2.x, p1.y}};
    }
    Vec2i shrink(Vec2i p) const {
        p.x = clamp(p.x, p1.x, p2.x);
        p.y = clamp(p.y, p1.y, p2.y);
        return p;
    }
    Vec2i center() const {
        Vec2i ret(p1 + p2);
        ret.x /= 2;
        ret.y /= 2;
        return ret;
    }
};

[[maybe_unused]] ostream &operator<<(ostream &out, const Rect &r) {
    return out << "(" << r.p1 << "), (" << r.p2 << ")";
}

template<class T, unsigned N>
struct SmallVec {
    struct RangeConstructTag {};

    unsigned size_ = 0;
    alignas(T) char data_[sizeof(T) * N];

    [[gnu::always_inline, gnu::pure, nodiscard]] T *data() noexcept {
        return reinterpret_cast<T*>(data_);
    }

    [[gnu::always_inline, gnu::pure, nodiscard]] const T *data() const noexcept {
        return reinterpret_cast<const T*>(data_);
    }

    [[gnu::always_inline, gnu::pure, nodiscard]] T &operator[](unsigned i) noexcept {
        return data()[i];
    }
    [[gnu::always_inline, gnu::pure, nodiscard]] const T &operator[](unsigned i) const noexcept {
        return data()[i];
    }

    [[gnu::always_inline, gnu::pure, nodiscard]] T &back() noexcept {
        return data()[size_ - 1];
    }
    [[gnu::always_inline, gnu::pure, nodiscard]] const T &back() const noexcept {
        return data()[size_ - 1];
    }

    [[gnu::always_inline, gnu::pure, nodiscard]] unsigned size() const {
        return size_;
    }

    [[gnu::always_inline, gnu::pure, nodiscard]] bool empty() const {
        return !size_;
    }

    [[gnu::always_inline, gnu::pure, nodiscard]] constexpr static unsigned capacity() {
        return N;
    }

    [[gnu::always_inline]] void resize(unsigned s) {
        if (s < size_) {
            for (unsigned i = s; i < size_; ++i) {
                data()[i].~T();
            }
        } else {
            for (unsigned i = size_; i < s; ++i) {
                new (data() + i) T{};
            }
        }
        size_ = s;
    }

    [[gnu::always_inline]] void resize(unsigned s, const T &val) {
        if (s < size_) {
            for (unsigned i = s; i < size_; ++i) {
                data()[i].~T();
            }
        } else {
            for (unsigned i = size_; i < s; ++i) {
                new (data() + i) T(val);
            }
        }
        size_ = s;
    }

    [[gnu::always_inline]] void clear() {
        for (unsigned i = 0; i < size_; ++i) {
            data()[i].~T();
        }
        size_ = 0;
    }

    template<class... Args>
    [[gnu::always_inline]] void emplace_back(Args&&... args) {
        new (&data()[size_]) T(std::forward<Args>(args)...);
        ++size_;
    }

    [[gnu::always_inline]] void push_back(T &&x) {
        emplace_back(std::move(x));
    }

    [[gnu::always_inline]] void push_back(const T &x) {
        emplace_back(x);
    }

    [[gnu::always_inline]] void pop_back() {
        --size_;
        data()[size_].~T();
    }

    SmallVec() = default;

    [[gnu::always_inline]] explicit SmallVec(unsigned n): size_(n) {
        for (unsigned i = 0; i < n; ++i)
            new (data() + i) T{};
    }

    [[gnu::always_inline]] explicit SmallVec(unsigned n, const T &val): size_(n) {
        for (unsigned i = 0; i < n; ++i)
            new (data() + i) T(val);
    }

    template<class Range, enable_if_t<is_integral_v<decltype(declval<Range>().size())>, bool> = true>
    [[gnu::always_inline]] SmallVec(Range &&other, RangeConstructTag = {}): size_(other.size()) {
        auto it = other.begin();
        for (unsigned i = 0; i < size_; ++i) {
            if constexpr (is_reference_v<Range> && !is_const_v<remove_reference_t<decltype(*it)>>) {
                new (data() + i) T(*it);
            } else {
                new (data() + i) T(move(*it));
            }
        }
    }

    [[gnu::always_inline]] SmallVec(initializer_list<T> l): SmallVec(l, RangeConstructTag{}) {}

    [[gnu::always_inline]] ~SmallVec() {
        for (unsigned i = 0; i < size_; ++i) {
            data()[i].~T();
        }
    }

    using iterator = T*;
    using const_iterator = const T*;

    [[gnu::always_inline, gnu::pure, nodiscard]] iterator begin() {
        return data();
    }

    [[gnu::always_inline, gnu::pure, nodiscard]] const_iterator begin() const {
        return data();
    }

    [[gnu::always_inline, gnu::pure, nodiscard]] iterator end() {
        return data() + size_;
    }

    [[gnu::always_inline, gnu::pure, nodiscard]] const_iterator end() const {
        return data() + size_;
    }

    [[gnu::always_inline]] void erase(iterator it) {
        while (it != prev(end())) {
            *it = std::move(*next(it));
            ++it;
        }
        it->~T();
        --size_;
    }
};

template<class T, T *buf, class Size = unsigned>
struct SmallPtr {
    Size pos = 0;

    SmallPtr() {}

    SmallPtr(nullptr_t): pos(0) {}

    SmallPtr(T *ptr): pos(ptr - buf + 1) {}

    explicit SmallPtr(Size p): pos(p + 1) {}

    T &operator*() const {
        return *(buf + pos - 1);
    }

    T *operator->() const {
        return buf + pos - 1;
    }

    operator bool() const {
        return pos != 0;
    }
};

// Constants

constexpr int COLORS = 4, TYPES = 3, SIZE = 10000, BASE_SCAN_RADIUS = 800, DRONE_SPEED = 600, LIGHT_SCAN_RADIUS = 2000, BATTAREY_FOR_LIGHT = 5;
constexpr int DRONE_EMERGENCY_SPEED = 300;
constexpr int NOISE_RADIUS = 1400, MAX_BATTAREY = 30, FISH_MAX_SPEED = 400;
constexpr int MONSTER_RADIUS = 500, MONSTER_SPEED = 540, MONSTER_EXTRA_RADIUS = 300;
constexpr int y_min[TYPES] = {2500, 5000, 7500};
constexpr int y_max[TYPES] = {5000, 7500, 10000};
constexpr int MAX_SCORE = 96;
constexpr int MAX_TICK = 200, SAVE_DEPTH = 500;

struct SmitsiMaxNode {
    static constexpr size_t MAX_NODES = 1'000'000;
    static constexpr size_t MOVES = COLORS * TYPES + 1;
    static SmitsiMaxNode nodes[MAX_NODES];
    static uint32_t allocated;

    using NodePtr = SmallPtr<SmitsiMaxNode, nodes, uint32_t>;

    static NodePtr newNode() {
        nodes[allocated] = {};
        ++allocated;
        return NodePtr(allocated - 1);
    }

    NodePtr parent = nullptr;
    array<int, MOVES> child_time{};
    array<NodePtr, MOVES> children{};
    double score = 0;
    int visits = 0;

    NodePtr create_child(int i) {
        children[i] = newNode();
        children[i]->parent = this;
        return children[i];
    }

    [[gnu::optimize("Os")]] void print(ostream &out, int t, int k = 0) {
        out << string(k, '\t') << "time = " << t << "score = " << score << " visits = " << visits << " (" << score / visits << ")\n";
        for (unsigned move = 0; move < SmitsiMaxNode::MOVES; ++move) {
            auto c = children[move];
            if (c) {
                out << string(k + 1, '\t') << move << ":\n";
                c->print(out, child_time[move], k + 1);
            }
        }
    }
};

SmitsiMaxNode SmitsiMaxNode::nodes[SmitsiMaxNode::MAX_NODES];
uint32_t SmitsiMaxNode::allocated = 0;

[[gnu::always_inline]] inline double uct(SmitsiMaxNode::NodePtr v, int total_visits) {
    return v->score / v->visits + sqrt(2 * log(total_visits) / v->visits);
}

template<class M, class V>
void gauss(M &mat, V &b) {
    size_t n = mat.size();
    size_t m = mat[0].size();
    for (size_t i = 0, l = 0; i < n; ++i, ++l) {
        size_t j = n;
        while (j == n && l < m) {
            j = find_if(mat.begin() + i, mat.end(), [l](auto &r) {
                return r[l] > 0.001;
            }) - mat.begin();
            l += (j == n);
        }
        if (l == m)
            break;
        if (j != i) {
            swap(mat[i], mat[j]);
            swap(b[i], b[j]);
        }
        b[i] /= mat[i][l];
        for (size_t j = l + 1; j < m; ++j) {
            mat[i][j] /= mat[i][l];
        }
        mat[i][l] = 1;
        for (size_t j = 0; j < n; ++j) {
            if (j == i)
                continue;
            if (mat[j][l] != 0) {
                auto mul = mat[j][l];
                for (size_t k = l; k < m; ++k) {
                    mat[j][k] = mat[j][k] - mat[i][k] * mul;
                }
                b[j] = b[j] - b[i] * mul;
            }
        }
    }
}

template<class M, class V1, class V2, class It, class MaxVals, class F>
pair<int, bool> find_les_solution(const M &a, const V1 &b, V2 &ans, It free_vars_first, It free_vars_last, const MaxVals &max_values, F &&check) {
    size_t n = a.size(), m = ans.size();
    if (free_vars_first == free_vars_last) {
        for (size_t i = 0, j = 0; i < n; ++i) {
            while (j < m && a[i][j] == 0)
                ++j;
            if (j == m)
                break;
            auto v = b[i];
            for (size_t k = j + 1; k < m; ++k) {
                v -= ans[k] * a[i][k];
            }
            if (v < 0 || v > max_values[j]) {
                return {0, 0};
            }
            ans[j] = v;
            ++j;
        }
        return {1, check(ans)};
    }
    int ret = 0;
    for (int i = 0; i <= max_values[*free_vars_first]; ++i) {
        ans[*free_vars_first] = i;
        auto [r, s] = find_les_solution(a, b, ans, next(free_vars_first), free_vars_last, max_values, check);
        ret += r;
        if (ret && !r || s)
            return {ret, s};
    }
    return {ret, 0};
}

REGISTER_CMD(Wait, bool);
REGISTER_CMD(Move, Vec2i, bool);

struct Fish {
    int id;
    int color, type;
};

struct Drone {
    int id;
    Vec2i pos;
    Vec2i vel;
    int emergency, battery;
};

struct MyDrone: public Drone {
    std::variant<CmdWait, CmdMove> cmd;
};

struct MovingFish {
    int id;
    Vec2i pos, vel;
    int time;
    bool estimated;
};

}

// struct LoggingBuffer : public streambuf {
//     ostream &log;
//     LoggingBuffer(ostream &log): log(log) {}
//     int last_char = traits_type::eof();
//     int underflow() override {
//         if (last_char == traits_type::eof()) {
//             last_char = cin.get();
//             log.put(last_char);
//             if (last_char == '\n')
//                 log.flush();
//         }
//         return last_char;
//     }
//     int uflow() override {
//         int ret = underflow();
//         last_char = traits_type::eof();
//         return ret;
//     }
// };

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
    // istream cin(new LoggingBuffer(cerr));
    mt19937 rnd(42);
    ios::sync_with_stdio(0);
    cin.tie(nullptr);
    int fishes_count;
    cin >> fishes_count;

    array<array<bool, TYPES>, COLORS> have{};
    array<array<bool, TYPES>, COLORS> foe_have{};
    set<int> have_id, foe_have_id;
    array<array<bool, TYPES>, COLORS> semi_have{};
    array<array<bool, TYPES>, COLORS> foe_semi_have{};
    set<int> semi_have_id, foe_semi_have_id;

    map<int, Fish> fishes;
    for (int i = 0; i < fishes_count; i++) {
        int id;
        cin >> id;
        fishes[id].id = id;
        cin >> fishes[id].color >> fishes[id].type;
    }

    int tick = 0;

    map<int, MovingFish> fish_cache;
    map<int, int> fish_check_time;
    map<int, Rect> fish_zone;

    set<int> saving;

    std::map<int, Drone> old_drones_data;

    // game loop
    while (1) {
        ++tick;

        int my_score;
        cin >> my_score;
        int foe_score;
        cin >> foe_score;
        int my_scan_count;
        cin >> my_scan_count;
        SmallVec<int, 12> my_scans(my_scan_count);
        for (auto &x : my_scans) {
            cin >> x;
        }
        int foe_scan_count;
        cin >> foe_scan_count;
        SmallVec<int, 12> foe_scans(foe_scan_count);
        for (auto &x : foe_scans) {
            cin >> x;
        }
        std::map<int, Drone&> drones;
        int my_drone_count;
        cin >> my_drone_count;
        SmallVec<MyDrone, 2> my_drones(my_drone_count);
        for (auto &d : my_drones) {
            auto &[id, pos, vel, e, b] = (Drone&) d;
            cin >> id >> pos >> e >> b;
            drones.emplace(id, d);
        }
        auto is_my = [&](int id) {
            for (auto &d : my_drones) {
                if (d.id == id) {
                    return 1;
                }
            }
            return 0;
        };
        int foe_drone_count;
        cin >> foe_drone_count;
        SmallVec<Drone, 2> foe_drones(foe_drone_count);
        for (auto &d : foe_drones) {
            auto &[id, pos, vel, e, b] = d;
            cin >> id >> pos >> e >> b;
            drones.emplace(id, d);
        }
        for (auto &[id, d] : drones) {
            if (old_drones_data.count(id)) {
                d.vel = d.pos - old_drones_data.at(id).pos;
            }
        }
        int drone_scan_count;
        cin >> drone_scan_count;
        struct DroneScan {
            int drone, fish;
        };
        SmallVec<DroneScan, 12 * 4> drone_scans(drone_scan_count);
        for (auto &[d, f] : drone_scans) {
            cin >> d >> f;
        }
        int visible_creature_count;
        cin >> visible_creature_count;
        for (int i = 0; i < visible_creature_count; i++) {
            MovingFish fish{};
            fish.time = tick;
            cin >> fish.id >> fish.pos >> fish.vel;
            fish_cache[fish.id] = fish;
        }

        int radar_blip_count;
        cin >> radar_blip_count;
        SmallVec<tuple<int, int, Vec2i>, 32 * 4> radare(radar_blip_count);
        for (auto &[d, f, dir] : radare) {
            string s;
            cin >> d >> f >> s;
            dir.y = (s[0] == 'T' ? -1 : +1);
            dir.x = (s[1] == 'L' ? -1 : +1);
        }

        // Logic

        for (auto &[id, f] : fish_cache) {
            if (f.time == tick)
                continue;
            int y1 = y_min[0];
            int y2 = SIZE;
            constexpr int x1 = 0;
            constexpr int x2 = SIZE;
            int t = fishes[f.id].type;
            if (t >= 0) {
                y1 = y_min[t];
                y2 = y_max[t];
            }
            if (f.pos.x + f.vel.x < x1 || f.pos.x + f.vel.x >= x2) {
                f.vel.x = -f.vel.x;
            }
            if (f.pos.y + f.vel.y < y1 || f.pos.y + f.vel.y >= y2) {
                f.vel.y = -f.vel.y;
            }
            f.pos += f.vel;
        }

        for (auto x : my_scans) {
            auto [_, c, t] = fishes[x];
            have[c][t] = 1;
            have_id.emplace(x);
        }

        for (auto x : foe_scans) {
            auto [_, c, t] = fishes[x];
            foe_have[c][t] = 1;
            foe_have_id.emplace(x);
        }

        semi_have_id = have_id;
        semi_have = have;

        for (auto [d, f] : drone_scans) {
            auto [_, c, t] = fishes[f];
            if (is_my(d)) {
                semi_have_id.emplace(f);
                semi_have[c][t] = 1;
            } else {
                bool is_new = foe_semi_have_id.emplace(f).second;
                if (is_new && fish_cache[f].time != tick) {
                    bool light = old_drones_data.count(d) && drones.at(d).battery < old_drones_data.at(d).battery;
                    fish_cache[f] = MovingFish{
                        .id = f,
                        .pos = Vec2i(drones.at(d).pos + drones.at(d).vel.normalize() * (light ? LIGHT_SCAN_RADIUS : BASE_SCAN_RADIUS)),
                        .vel = {},
                        .time = tick,
                        .estimated = 1,
                    };
                }
                foe_semi_have[c][t] = 1;
            }
        }

        // Expand zones
        for (auto &[id, fish] : fishes) {
            int t = fish.type;
            int y1 = (t == -1 ? y_min[0] : y_min[t]);
            int y2 = (t == -1 ? SIZE : y_max[t]);
            Rect max_zone{{0, y1}, {SIZE, y2}};
            if (!fish_zone.count(id)) {
                fish_zone[id] = max_zone;
            } else {
                auto &zone = fish_zone[id];
                zone.p1.x -= FISH_MAX_SPEED;
                zone.p1.y -= FISH_MAX_SPEED;
                zone.p2.x += FISH_MAX_SPEED;
                zone.p2.y += FISH_MAX_SPEED;
                zone = zone.intersect(max_zone);
            }
        }

        for (auto [d, f, dir] : radare) {
            Rect bounds;
            if (dir.x < 0) {
                bounds.p1.x = 0;
                bounds.p2.x = drones.at(d).pos.x;
            } else {
                bounds.p1.x = drones.at(d).pos.x;
                bounds.p2.x = SIZE;
            }
            if (dir.y < 0) {
                bounds.p1.y = 0;
                bounds.p2.y = drones.at(d).pos.y;
            } else {
                bounds.p1.y = drones.at(d).pos.y;
                bounds.p2.y = SIZE;
            }
            fish_zone[f] = fish_zone[f].intersect(bounds);
        }

        auto is_visble_by = [&](const MovingFish &fish, const Drone &drone) {
            bool light = (old_drones_data.count(drone.id) && old_drones_data.at(drone.id).battery > drone.battery);
            int radius = (light ? LIGHT_SCAN_RADIUS : BASE_SCAN_RADIUS);
            if (fishes[fish.id].type < 0)
                radius += MONSTER_EXTRA_RADIUS;
            return (fish.pos - drone.pos).len() < radius;
        };

        for (auto &d : my_drones) {
            for (auto &[id, fish] : fish_cache) {
                if (is_visble_by(fish, d)) {
                    fish_check_time[id] = tick;
                }
            }
        }

        // Smitsimax

        array<Fish, COLORS * TYPES> flat_fishes;
        array<Vec2i, COLORS * TYPES> flat_pos;
        array<bool, COLORS * TYPES> invalid{};
        {
            auto it_f = flat_fishes.begin();
            auto it_p = flat_pos.begin();
            auto it_i = invalid.begin();
            for (auto &[id, fish] : fishes) {
                if (fish.type >= 0) {
                    *it_f++ = fish;
                    *it_p++ = fish_zone[id].center();
                    *it_i++ = find_if(radare.begin(), radare.end(), [&](auto &r) {
                        return get<1>(r) == id;
                    }) == radare.end();
                }
            }
        }

        array<int, 4> drone_to_id;
        array<Vec2i, 4> initial_coords;

        {
            auto it_id = drone_to_id.begin();
            auto it_coord = initial_coords.begin();
            for (unsigned i = 0; i < my_drones.size(); ++i) {
                *it_id++ = my_drones[i].id;
                *it_id++ = foe_drones[i].id;
                *it_coord++ = my_drones[i].pos;
                *it_coord++ = foe_drones[i].pos;
            }
        }

        array<array<int, TYPES>, 2> intitial_found_of_type{};
        array<array<int, COLORS>, 2> intitial_found_of_color{};
        array<SmallVec<int, COLORS * TYPES>, 4> initial_drone_scans{};
        array<array<bool, COLORS>, 2> initial_color_bonus_possible{};
        array<array<bool, TYPES>, 2> initial_type_bonus_possible{};
        array<array<bool, COLORS * TYPES>, 2> initial_fish_bonus_possible{};

        for (auto [d, f] : drone_scans) {
            int drone_index = find(drone_to_id.begin(), drone_to_id.end(), d) - drone_to_id.begin();
            int fish_index = find_if(flat_fishes.begin(), flat_fishes.end(), [f](auto &fish) {
                return fish.id == f;
            }) - flat_fishes.begin();
            initial_drone_scans[drone_index].emplace_back(fish_index);
        }

        for (auto x : my_scans) {
            auto [_, c, t] = fishes[x];
            ++intitial_found_of_color[0][c];
            ++intitial_found_of_type[0][t];
        }

        for (auto x : foe_scans) {
            auto [_, c, t] = fishes[x];
            ++intitial_found_of_color[1][c];
            ++intitial_found_of_type[1][t];
        }

        array<SmitsiMaxNode::NodePtr, 4> roots;

        SmitsiMaxNode::allocated = 0;
        for (auto &p : roots) {
            p = SmitsiMaxNode::newNode();
        }

        array<array<bool, COLORS * TYPES>, 2> initial_useless_fish{};
        array<int, 2> initial_possible_score = {my_score, foe_score};

        {
            array<array<int, COLORS>, 2> color_cnt{};
            array<array<int, TYPES>, 2> type_cnt{};
            for (int i = 0; i < COLORS * TYPES; ++i) {
                auto [_, c, t] = flat_fishes[i];
                initial_fish_bonus_possible[0][i] = (!invalid[i] || semi_have[c][t]);
                initial_fish_bonus_possible[1][i] = (!invalid[i] || foe_semi_have[c][t]);
                if (!have[c][t] && initial_fish_bonus_possible[0][i])
                    initial_possible_score[0] += t + 1;
                if (!foe_have[c][t] && initial_fish_bonus_possible[1][i])
                    initial_possible_score[1] += t + 1;
                for (int p = 0; p < 2; ++p) {
                    color_cnt[p][c] += (initial_fish_bonus_possible[p][i]);
                    type_cnt[p][t] += (initial_fish_bonus_possible[p][i]);
                    initial_fish_bonus_possible[p][i] &= !have[c][t] && !foe_have[c][t];
                    initial_possible_score[p] += (t + 1) * initial_fish_bonus_possible[p][i];
                }
                if (semi_have[c][t] || invalid[i])
                    initial_useless_fish[0][i] = 1;
                if (foe_semi_have[c][t] || invalid[i])
                    initial_useless_fish[1][i] = 1;
            }
            cerr << "possible_score: " << initial_possible_score[0] << " " << initial_possible_score[1] << "\n";
            cerr << "color_cnt:\n";
            for (auto &v : color_cnt) {
                for (auto x : v) {
                    cerr << x << " ";
                }
                cerr << "\n";
            }
            cerr << "type_cnt:\n";
            for (auto &v : type_cnt) {
                for (auto x : v) {
                    cerr << x << " ";
                }
                cerr << "\n";
            }
            for (int p = 0; p < 2; ++p) {
                for (int i = 0; i < COLORS; ++i) {
                    if (color_cnt[p][i] == TYPES) {
                        initial_color_bonus_possible[p][i] = max(intitial_found_of_color[0][i], intitial_found_of_color[1][i]) < TYPES;
                        initial_possible_score[p] += TYPES * ((intitial_found_of_color[p][i] < TYPES) + initial_color_bonus_possible[p][i]);
                    }
                }
                for (int i = 0; i < TYPES; ++i) {
                    if (type_cnt[p][i] == COLORS) {
                        initial_type_bonus_possible[p][i] = max(intitial_found_of_type[0][i], intitial_found_of_type[1][i]) < COLORS;
                        initial_possible_score[p] += COLORS * ((intitial_found_of_type[p][i] < COLORS) + initial_type_bonus_possible[p][i]);
                    }
                }
            }
        }

        array<int, 4> initial_move_time;

        for (int i = 0; i < 4; ++i) {
            auto &d = drones.at(drone_to_id[i]);
            initial_move_time[i] = tick * DRONE_SPEED;
            if (d.emergency) {
                initial_coords[i].y = 0;
                initial_move_time[i] += (d.pos.y + DRONE_EMERGENCY_SPEED - 1) / DRONE_EMERGENCY_SPEED * DRONE_SPEED;
            }
        }

        auto start = clock();
        #ifndef DETERMINISTIC_MODE
        constexpr long TIME_LIMIT = CLOCKS_PER_SEC * 47 / 1000;
        #endif
        array<double, 2> sum_score{};
        array<int, 2> min_score = {200, 200}, max_score = {0, 0};
        int iterations = 0;
        #ifdef DETERMINISTIC_MODE
        while (iterations < 7000) {
        #else
        while (clock() - start < TIME_LIMIT) {
        #endif
            for (int i = 0; i < 100; ++i) {
                auto cur = roots;
                auto used = initial_useless_fish;

                array<int, 4> last_move;
                for (int i = 0; i < 4; ++i) {
                    if (initial_drone_scans[i].empty()) {
                        last_move[i] = 0;
                    } else {
                        last_move[i] = -1;
                    }
                }
                array<array<bool, COLORS * TYPES>, 2> saved = {};
                for (int i = 0; i < COLORS * TYPES; ++i) {
                    auto [_, c, t] = flat_fishes[i];
                    if (have[c][t])
                        saved[0][i] = 1;
                    if (foe_have[c][t])
                        saved[1][i] = 1;
                }
                array<int, 2> score = {my_score, foe_score};
                auto found_of_type = intitial_found_of_type;
                auto found_of_color = intitial_found_of_color;
                auto drone_scans = initial_drone_scans;
                auto coords = initial_coords;

                auto type_bonus_possible = initial_type_bonus_possible;
                auto color_bonus_possible = initial_color_bonus_possible;
                auto fish_bonus_possible = initial_fish_bonus_possible;
                auto possible_score = initial_possible_score;

                auto move_time = initial_move_time;

                while (1) {
                    auto drone = min_element(move_time.begin(), move_time.end()) - move_time.begin();
                    if (move_time[drone] >= MAX_TICK * DRONE_SPEED)
                        break;

                    int player = drone & 1;

                    if (last_move[drone] == 0) {
                        for (auto x : drone_scans[drone]) {
                            if (!saved[player][x]) {
                                saved[player][x] = 1;
                                auto &[_, c, t] = flat_fishes[x];
                                ++found_of_color[player][c];
                                ++found_of_type[player][t];
                                possible_score[!player] -= fish_bonus_possible[!player][x] * (t + 1);
                                fish_bonus_possible[!player][x] = 0;
                                int add = t + 1;
                                if (!saved[!player][x]) {
                                    add *= 2;
                                }
                                bool foe_have_c = (found_of_color[!player][c] == TYPES);
                                bool foe_have_t = (found_of_type[!player][t] == COLORS);
                                if (found_of_color[player][c] == TYPES) {
                                    add += TYPES * (foe_have_c ? 1 : 2);
                                    possible_score[!player] -= color_bonus_possible[!player][c] * TYPES;
                                    color_bonus_possible[!player][c] = 0;
                                }
                                if (found_of_type[player][t] == COLORS) {
                                    add += COLORS * (foe_have_t ? 1 : 2);
                                    possible_score[!player] -= type_bonus_possible[!player][t] * COLORS;
                                    type_bonus_possible[!player][t] = 0;
                                }
                                score[player] += add;// * exp((tick - move_time[i]) * 0.007);
                                // if (iterations == 0 && i == 0 && player == 1) {
                                //     cerr << c << " " << t << " +" << add << "\n";
                                // }
                            }
                        }
                        if (score[player] > possible_score[!player]) {
                            break;
                        }
                        drone_scans[drone].clear();
                    }

                    SmitsiMaxNode::NodePtr best = nullptr;
                    int best_move = 0;

                    if (cur[drone]->visits == 0) {
                        for (unsigned i = (last_move[drone] == 0); i < SmitsiMaxNode::MOVES; ++i) {
                            if (i && used[player][i - 1])
                                continue;
                            if (i == 0) {
                                cur[drone]->child_time[i] = move_time[drone] + coords[drone].y - SAVE_DEPTH;
                            } else {
                                cur[drone]->child_time[i] = move_time[drone] + int((coords[drone] - flat_pos[i - 1]).len());
                            }
                        }
                    }

                    if (!best) {
                        for (unsigned i = 0; i < SmitsiMaxNode::MOVES; ++i) {
                            auto ch = cur[drone]->children[i];
                            if (ch && ch->score == ch->visits) {
                                best_move = i;
                                best = ch;
                                break;
                            }
                        }
                    }

                    if (!best && last_move[drone] != 0 && !cur[drone]->children[0]) {
                        best = cur[drone]->create_child(best_move = 0);
                    }

                    // int max_possible_time = MAX_TICK * DRONE_SPEED;
                    // array<bool, COLORS * TYPES> my_fishes{};
                    // for (auto x : drone_scans[drone])
                    //     my_fishes[x] = 1;

                    // int my_min_time = coords[drone].y - SAVE_DEPTH;

                    // #pragma GCC unroll 2
                    // for (int foe_drone = !player; foe_drone < 4; foe_drone += 2) {
                    //     bool intersect = 0;
                    //     for (auto x : drone_scans[foe_drone]) {
                    //         if (my_fishes[x]) {
                    //             intersect = 1;
                    //             break;
                    //         }
                    //     }
                    //     if (intersect) {
                    //         int foe_time = move_time[foe_drone];
                    //         if (last_move[foe_drone] > 0) {
                    //             foe_time += (coords[foe_drone].y - SAVE_DEPTH);
                    //         }
                    //         if (foe_time > my_min_time) {
                    //             max_possible_time = min(max_possible_time, foe_time);
                    //         }
                    //     }
                    // }

                    if (!best) {
                        int best_fish = -1;
                        int best_dist = SIZE * 2;
                        int best_type = 0;
                        for (unsigned i = 0; i < COLORS * TYPES; ++i) {
                            if (used[player][i] || cur[drone]->child_time[i + 1] >= MAX_TICK * DRONE_SPEED || cur[drone]->children[i + 1])
                                continue;
                            int dist = (flat_pos[i] - coords[drone]).len();
                            int type = flat_fishes[i].type;
                            if (pair(type, -dist) > pair(best_type, -best_dist)) {
                                best_dist = dist;
                                best_type = type;
                                best_fish = i;
                            }
                        }

                        if (best_fish != -1) {
                            best = cur[drone]->create_child(best_move = best_fish + 1);
                        } else if (!cur[drone]->children[0] && last_move[drone] != 0) {
                            best = cur[drone]->create_child(best_move = 0);
                        }
                    }
                    
                    if (!best) {
                        bool unused = 0;
                        double best_uct = 0;
                        for (size_t i = (last_move[drone] == 0); i < SmitsiMaxNode::MOVES; ++i) {
                            if (i && saved[player][i - 1])
                                continue;
                            bool i_unused = (i == 0 || !used[player][i - 1]);
                            double u = uct(cur[drone]->children[i], cur[drone]->visits);
                            if (!best || pair(i_unused, u) > pair(unused, best_uct)) {
                                best_uct = u;
                                unused = i_unused;
                                best = cur[drone]->children[best_move = i];
                            }
                        }
                    }

                    if (!best) {
                        move_time[drone] = MAX_TICK * DRONE_SPEED;
                        continue;
                    }

                    cur[drone] = best;
                    last_move[drone] = best_move;
                    if (best_move == 0) {
                        move_time[drone] += max(DRONE_SPEED, coords[drone].y - SAVE_DEPTH);
                        coords[drone].y = SAVE_DEPTH;
                    } else {
                        Vec2d d = (coords[drone] - flat_pos[best_move - 1]);
                        if (d.len() > LIGHT_SCAN_RADIUS) {
                            d = d.normalize() * LIGHT_SCAN_RADIUS;
                        }
                        Vec2i new_coord = flat_pos[best_move - 1] + Vec2i(d);
                        move_time[drone] += max<int>((new_coord - coords[drone]).len(), 1);
                        coords[drone] = new_coord;
                        drone_scans[drone].emplace_back(best_move - 1);
                        used[player][best_move - 1] = 1;
                    }
                }

                array<double, 2> normalized_scores;
                for (int i = 0; i < 2; ++i) {
                    if (score[i] > possible_score[!i])
                        normalized_scores[i] = 1;
                    else
                        normalized_scores[i] = 0.5 * possible_score[i] / (possible_score[0] + possible_score[1]);
                }

                for (int i = 0; i < 2; ++i) {
                    sum_score[i] += score[i];
                    min_score[i] = min(min_score[i], score[i]);
                    max_score[i] = max(max_score[i], score[i]);
                }

                for (int drone = 0; drone < 4; ++drone) {
                    int player = drone & 1;
                    while (cur[drone]) {
                        // cerr << drone << ": " << cur[drone].pos << "\n";
                        cur[drone]->visits++;
                        cur[drone]->score += normalized_scores[player];
                        cur[drone] = cur[drone]->parent;
                    }
                }
            }
            iterations += 100;
        }
        auto end = clock();

        cerr << "possible_score: " << initial_possible_score[0] << " " << initial_possible_score[1] << "\n";
        cerr << "avg_score: " << sum_score[0] / iterations << " " << sum_score[1] / iterations << "\n";
        cerr << "min_score: " << min_score[0] << " " << min_score[1] << "\n";
        cerr << "max_score: " << max_score[0] << " " << max_score[1] << "\n";
        cerr << "Iterations: " << iterations << " / " << static_cast<double>(end - start) * 1000.0 / CLOCKS_PER_SEC << " ms\n";
        for (int i = 0; i < 2; ++i) {
            assert(max_score[i] <= initial_possible_score[i]);
        }
        cerr << "Nodes allocated: " << SmitsiMaxNode::allocated << "\n";
        // if (SmitsiMaxNode::allocated < 50) {
        //     int i = 0;
        //     for (auto &[_, c, t] : flat_fishes) {
        //         ++i;
        //         cerr << i << ": " << c << " " << t << "\n";
        //     }
        //     for (int i = 0; i < 4; ++i)
        //         roots[i]->print(cerr, initial_move_time[i]);
        // }

        for (int i = 0; i < 4; ++i) {
            cerr << "Drone " << drone_to_id[i] << " root score: " << roots[i]->score / roots[i]->visits << "\n";
        }

        for (int i = 0; i < 4; i += 2) {
            auto &cur_drone = my_drones[i / 2];
            double best_score = -1;
            int best_move = 0;
            int best_visits = 0;
            for (size_t m = 0; m < SmitsiMaxNode::MOVES; ++m) {
                auto ch = roots[i]->children[m];
                if (ch && ch->visits && ch->score / ch->visits > best_score) {
                    best_score = ch->score / ch->visits;
                    best_move = m;
                    best_visits = ch->visits;
                }
            }
            if (best_move == 0) {
                cerr << "Best move for drone " << cur_drone.id << ": go up, score: " << best_score << " (" << best_visits << ")\n";
                auto target = cur_drone.pos;
                target.y = SAVE_DEPTH - 1;
                cur_drone.cmd = CmdMove{{
                    target,
                    0,
                }};
            } else {
                --best_move;
                auto [id, c, t] = flat_fishes[best_move];
                cerr << "Best move for drone " << cur_drone.id << ": go to fish " << id << "(c: " << c << ", t: " << t << ") at (" << flat_pos[best_move];
                cerr << ") score: " << best_score << " (" << best_visits << ")\n";
                double dist = (fish_zone[id].center() - cur_drone.pos).len();
                cur_drone.cmd = CmdMove{{
                    flat_pos[best_move],
                    dist > BASE_SCAN_RADIUS && dist < LIGHT_SCAN_RADIUS && cur_drone.battery >= BATTAREY_FOR_LIGHT || cur_drone.battery == MAX_BATTAREY && cur_drone.pos.y >= y_min[0],
                }};
            }
        }

        // set<int> scared_fishes;

        // for (int i = 0; i < 2; ++i) {
        //     auto &cur_drone = my_drones[i];
        //     [[maybe_unused]] auto &other_my_drone = my_drones[i ^ 1];

        //     // Scare fishes

        //     auto is_real_and_visble = [&](const MovingFish &fish) {
        //         return !fish.estimated && fish.time == tick && is_visble_by(fish, cur_drone);
        //     };

        //     int min_dist = SIZE;
        //     int best_fish = -1;
        //     Vec2i scare_dir{};

        //     for (auto &[_, fish] : fish_cache) {
        //         if (fishes[fish.id].type < 0)
        //             continue;
        //         if (!is_real_and_visble(fish))
        //             continue;
        //         if (foe_semi_have_id.count(fish.id) || scared_fishes.count(fish.id))
        //             continue;
        //         if (min(fish.pos.x, SIZE - fish.pos.x) > min_dist)
        //             continue;
        //         cerr << "Drone " << cur_drone.id << " trys to scare fish " << fish.id << "\n";
        //         Vec2i dirs[] = {{-1, 0}, {1, 0}};
        //         if (fish.pos.x > SIZE / 2) {
        //             swap(dirs[0], dirs[1]);
        //         }
        //         for (auto dir : dirs) {
        //             bool ok = 1;
        //             for (auto &d : foe_drones) {
        //                 if ((d.pos - fish.pos) * dir > 0) {
        //                     ok = 0;
        //                     break;
        //                 }
        //             }
        //             if (ok) {
        //                 best_fish = fish.id;
        //                 min_dist = min(fish.pos.x, SIZE - fish.pos.x);
        //                 scare_dir = dir;
        //                 break;
        //             }
        //         }
        //     };

        //     if (best_fish != -1) {
        //         scared_fishes.emplace(best_fish);
        //         bool light = (fish_cache[best_fish].pos - cur_drone.pos).len() >= BASE_SCAN_RADIUS;
        //         cur_drone.cmd = CmdMove{{fish_cache[best_fish].pos - scare_dir * (BASE_SCAN_RADIUS - FISH_MAX_SPEED - 10), light}};
        //         break;
        //     }

        //     // Check save needed

        //     bool need_to_save = 0, nothing_to_save = 1;

        //     for (auto [d, c] : drone_scans) {
        //         if (d != cur_drone.id)
        //             continue;
        //         if (have_id.count(c) || foe_have_id.count(c))
        //             continue;
        //         nothing_to_save = 0;
        //         for (auto [d, f] : drone_scans) {
        //             if (is_my(d) || f != c)
        //                 continue;
        //             if (drones.at(d).pos.y <= cur_drone.pos.y + DRONE_SPEED * 2) {
        //                 need_to_save = 1;
        //                 break;
        //             }
        //         }
        //         if (need_to_save)
        //             break;
        //         // cerr << "Drone " << cur_drone.id << " fish " << c << ": save not needed\n";
        //     }

        //     if (need_to_save) {
        //         saving.emplace(cur_drone.id);
        //     } else if (nothing_to_save) {
        //         saving.erase(cur_drone.id);
        //     }

        //     if (saving.count(cur_drone.id)) {
        //         auto target = cur_drone.pos;
        //         target.y = 0;
        //         cur_drone.cmd = CmdMove{{target, 0}};
        //         continue;
        //     }

        //     // Regular scan
        //     // Vec2i target;
        //     // {
        //     //     constexpr auto W = (SIZE / DRONE_SPEED);
        //     //     constexpr auto H = (SIZE / (BASE_SCAN_RADIUS * 2));
        //     //     target.x = BASE_SCAN_RADIUS + (tick / W % H) * BASE_SCAN_RADIUS * 2;
        //     //     target.y = DRONE_SPEED / 2 + (tick % W) * DRONE_SPEED;
        //     //     if (tick / (H * W) % 2 == i) {
        //     //         target.x = SIZE - target.x;
        //     //     }
        //     //     if (tick / W % H % 2 == i) {
        //     //         target.y = SIZE - target.y;
        //     //     }
        //     // }

        //     // double scan_radius;

        //     // if (cur_drone.battery >= BATTAREY_FOR_LIGHT) {
        //     //     scan_radius = LIGHT_SCAN_RADIUS;
        //     // } else {
        //     //     scan_radius = BASE_SCAN_RADIUS;
        //     // }

        //     // double best_score = 0;
        //     // bool light = 0;

        //     // for (auto &[id, fish] : fish_cache) {
        //     //     if (fish_check_time[id] >= fish.time)
        //     //         continue;
        //     //     auto [_, c, t] = fishes[id];
        //     //     if (semi_have_id.count(id))
        //     //         continue;
        //     //     double score = (t + 1) * 2;
        //     //     if (!foe_have_id.count(id)) {
        //     //         score *= 2;
        //     //     }
        //     //     int c_cnt = 0, t_cnt = 0;
        //     //     bool foe_have_c = 1, foe_have_t = 1;
        //     //     for (int i = 0; i < COLORS; ++i) {
        //     //         t_cnt += !semi_have[i][t];
        //     //         foe_have_t &= foe_have[i][t];
        //     //     }
        //     //     for (int i = 0; i < TYPES; ++i) {
        //     //         c_cnt += !semi_have[c][i];
        //     //         foe_have_c &= foe_have[c][i];
        //     //     }
        //     //     if (c_cnt == 1) {
        //     //         score += 3 * (foe_have_c ? 1 : 2);
        //     //     }
        //     //     if (t_cnt == 1) {
        //     //         score += 4 * (foe_have_t ? 1 : 2);
        //     //     }
        //     //     double dist = (cur_drone.pos - fish.pos).len();
        //     //     if (dist < scan_radius) {
        //     //         if (dist > BASE_SCAN_RADIUS)
        //     //             light = 1;
        //     //         continue;
        //     //     }
        //     //     int time = ceil((dist - scan_radius) / DRONE_SPEED);
        //     //     score /= time;
        //     //     if (score > best_score) {
        //     //         best_score = score;
        //     //         target = fish.pos;
        //     //     }
        //     // }

        //     // cerr << "Drone " << cur_drone.id << ": best_score = " << best_score << "\n";

        //     cur_drone.cmd = CmdMove{{
        //         target, 0
        //         //cur_drone.battery == MAX_BATTAREY || rects[best_rect[i]].contains(cur_drone.pos),
        //     }};
        // }

        // Handle monsters
        for (auto &cur_drone : my_drones) {
            constexpr auto DANGER_RADIUS = MONSTER_RADIUS;
            if (auto cmd = get_if<CmdMove>(&cur_drone.cmd)) {
                auto &[target, light] = *static_cast<CmdMove::Base*>(cmd);
                SmallVec<MovingFish, 10> monsters;
                for (auto [id, m] : fish_cache) {
                    if (fishes[id].type < 0 && fish_check_time[id] <= m.time) {
                        monsters.emplace_back(m);
                    }
                }
                sort(monsters.begin(), monsters.end(), [&](auto &a, auto &b) {
                    return (a.pos - cur_drone.pos).lenSq() > (b.pos - cur_drone.pos).lenSq();
                });
                if (cur_drone.id < 2) {
                    for (auto &m : monsters) {
                        cerr << "Monster " << m.id << " at " << m.pos << " time = " << m.time << " checked at " << fish_check_time[m.id] << "\n";
                    }
                }
                for (auto &m : monsters) {
                    Vec2d s = m.vel;
                    if ((m.pos - cur_drone.pos).len() <= BASE_SCAN_RADIUS) {
                        s = (cur_drone.pos - m.pos).normalize() * MONSTER_SPEED;
                    }
                    Vec2i new_m_pos(m.pos + s);
                    if ((new_m_pos.x - cur_drone.pos.x) * (m.pos.x - cur_drone.pos.x) < 0) {
                        new_m_pos.x = cur_drone.pos.x;
                    }
                    if ((new_m_pos.y - cur_drone.pos.y) * (m.pos.y - cur_drone.pos.y) < 0) {
                        new_m_pos.y = cur_drone.pos.y;
                    }
                    cerr << "new pos = " << new_m_pos << " dist: " << (new_m_pos - cur_drone.pos).len() << "\n";
                    // if ((new_m_pos - cur_drone.pos).len() < LIGHT_SCAN_RADIUS && (new_m_pos - cur_drone.pos).len() > BASE_SCAN_RADIUS) {
                    //     light = 0;
                    // }
                    if ((cur_drone.pos - new_m_pos).len() <= MONSTER_RADIUS) {
                        if (new_m_pos == cur_drone.pos) {
                            new_m_pos = m.pos;
                        }
                        target = cur_drone.pos + Vec2i((cur_drone.pos - new_m_pos).normalize() * DRONE_SPEED);
                        if (target.x < 0 || target.x >= SIZE || target.y < 0 || target.y >= SIZE) {
                            Vec2d d = target - cur_drone.pos;
                            d.x += copysign(DRONE_SPEED, d.x);
                            d.y += copysign(DRONE_SPEED, d.y);
                            target = Vec2i(cur_drone.pos + d);
                        }
                        cerr << "Monster too close, go to " << target << "\n";
                        break;
                    }
                    // if ((cur_drone.pos - new_m_pos).len() <= DANGER_RADIUS) {
                    //     Vec2d d = cur_drone.pos - new_m_pos;
                    //     d /= d.len();
                    //     cerr << "Drone " << cur_drone.id << ": (" << target << ") -> ";
                    //     target = Vec2i(cur_drone.pos + d * DRONE_SPEED);
                    //     if (target.x < 0 || target.x >= SIZE || target.y < 0 || target.y >= SIZE) {
                    //         d.x += copysign(DRONE_SPEED, d.x);
                    //         d.y += copysign(DRONE_SPEED, d.y);
                    //         target = Vec2i(cur_drone.pos + d * DRONE_SPEED);
                    //     }
                    //     cerr << "(" << target << ") monster: " << m.id << "\n";
                    //     continue;
                    // }
                    if ((target - new_m_pos).len() <= DANGER_RADIUS) {
                        Vec2d d = target - new_m_pos;
                        d /= d.len();
                        cerr << "Drone " << cur_drone.id << ": (" << target << ") -> ";
                        target = Vec2i(new_m_pos + d * DANGER_RADIUS);
                        // if (target.x < 0 || target.x >= SIZE || target.y < 0 || target.y >= SIZE) {
                        //     d = target - cur_drone.pos;
                        //     d.x += copysign(DRONE_SPEED, d.x);
                        //     d.y += copysign(DRONE_SPEED, d.y);
                        //     target = Vec2i(cur_drone.pos + d);
                        // }
                        cerr << "(" << target << ") monster: " << m.id << " case 1\n";
                        // continue;
                    }
                    auto target_dir = target - cur_drone.pos;
                    auto monster_dir = new_m_pos - cur_drone.pos;
                    cerr << "monster_dir = " << monster_dir << "\n";
                    cerr << "target_dir = " << target_dir << "\n";
                    if (monster_dir * target_dir >= target_dir * target_dir || monster_dir * target_dir <= 0)
                        continue;
                    auto dist = abs(cross(monster_dir, target_dir) / target_dir.len());
                    if (dist > DANGER_RADIUS)
                        continue;
                    double free_zone_radius = (cur_drone.pos - new_m_pos).len() < BASE_SCAN_RADIUS ? MONSTER_RADIUS : BASE_SCAN_RADIUS;
                    double sn = (free_zone_radius + 10) / (cur_drone.pos - new_m_pos).len();
                    double cs = sqrt(1 - sn * sn);

                    Vec2d rot(monster_dir.x * cs - monster_dir.y * sn, monster_dir.x * sn + monster_dir.y * cs);
                    if (signbit(cross(rot, monster_dir)) != signbit(cross(target_dir, monster_dir))) {
                        rot = Vec2d(monster_dir.x * cs + monster_dir.y * sn, monster_dir.y * cs - monster_dir.x * sn);
                    }
                    cerr << "Drone " << cur_drone.id << ": (" << target << ") -> ";
                    target = cur_drone.pos + Vec2i(rot);
                    if (target.x < 0 || target.x >= SIZE || target.y < 0 || target.y >= SIZE) {
                        Vec2d d = target - cur_drone.pos;
                        d.x += copysign(DRONE_SPEED, d.x);
                        d.y += copysign(DRONE_SPEED, d.y);
                        target = Vec2i(cur_drone.pos + d);
                    }
                    cerr << "(" << target << ") monster: " << m.id << " case 2\n";
                }
            }
        }

        // Shrink commands

        for (auto &d : my_drones) {
            if (auto cmd = get_if<CmdMove>(&d.cmd)) {
                auto &[target, light] = *static_cast<CmdMove::Base*>(cmd);
                auto s = (target - d.pos);
                if (s.len() > DRONE_SPEED) {
                    Vec2i new_target(d.pos + s.normalize() * DRONE_SPEED);
                    if (new_target.x >= 0 && new_target.x < SIZE && new_target.y >= 0 && new_target.y <= SIZE)
                        continue;
                }
                target.x = clamp(target.x, 0, SIZE - 1);
                target.y = clamp(target.y, 0, SIZE - 1);
            }
        }

        // Update game data

        old_drones_data = decltype(old_drones_data)(drones.begin(), drones.end());

        // Apply commands
        for (auto &d : my_drones) {
            visit([](auto cmd) {
                for (auto x : decltype(cmd)::mnemonic) {
                    cout << (char) toupper(x);
                }
                [&]<size_t... i>(index_sequence<i...>) {
                    ((cout << " " << get<i>(cmd)), ...);
                }(make_index_sequence<tuple_size_v<typename decltype(cmd)::Base>>{});
                cout << endl;
            }, d.cmd);
        }
    }
}
