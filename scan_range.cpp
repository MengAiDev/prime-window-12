#include <bits/stdc++.h>
#include <atomic>
#include <thread>
using namespace std;
typedef unsigned long long u64;
typedef unsigned int u32;

static u64 L    = 1011128158584751ULL;
static u64 R    = 10098097238186292ULL;
static const u64 WIN  = 2003ULL;              // interval [N, N+WIN]
static const u64 SMAX = 100500000ULL;         // sqrt(R) + margin

static vector<u64> base_primes;

// ---------------- base primes up to SMAX ----------------
void gen_base_primes() {
    u64 n_odds = (SMAX - 3) / 2 + 1;
    u64 nw = (n_odds + 63) / 64;
    vector<u64> base(nw, 0ULL);
    auto setc = [&](u64 x) {
        u64 idx = (x - 3) / 2;
        base[idx >> 6] |= 1ULL << (idx & 63);
    };
    auto isc = [&](u64 x) -> bool {
        u64 idx = (x - 3) / 2;
        return ((base[idx >> 6] >> (idx & 63)) & 1ULL) != 0;
    };
    for (u64 p = 3; p * p <= SMAX; p += 2)
        if (!isc(p))
            for (u64 m = p * p; m <= SMAX; m += 2 * p)
                setc(m);
    base_primes.reserve(6000000);
    base_primes.push_back(2);
    for (u64 p = 3; p <= SMAX; p += 2)
        if (!isc(p)) base_primes.push_back(p);
}

// ---------------- deterministic Miller-Rabin ----------------
static inline u64 mulmod(u64 a, u64 b, u64 mod) { return (u64)((__uint128_t)a * b % mod); }
static inline u64 powmod(u64 a, u64 d, u64 mod) {
    u64 r = 1;
    while (d) { if (d & 1) r = mulmod(r, a, mod); a = mulmod(a, a, mod); d >>= 1; }
    return r;
}
bool is_prime_mr(u64 n) {
    if (n < 2) return false;
    for (u64 p : {2ULL,3ULL,5ULL,7ULL,11ULL,13ULL,17ULL,19ULL,23ULL,29ULL,31ULL,37ULL})
        if (n % p == 0) return n == p;
    u64 d = n - 1, s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (u64 a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        a %= n; if (a == 0) continue;
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool ok = false;
        for (u64 r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { ok = true; break; } }
        if (!ok) return false;
    }
    return true;
}

// ---------------- incremental segmented sieve ----------------
struct SegSieve {
    u64 BS;
    u64 lo_odd, hi_odd;
    u32 n_odd;
    vector<u64> d;          // per-prime offset of next odd multiple to mark (even, >= 0)
    vector<u32> bsp;        // BS % (2*p)
    vector<unsigned char> comp;
    vector<u64> pr;

    void set_block(u64 S, u64 maxN) {
        lo_odd = (S & 1) ? S : S + 1;
        hi_odd = maxN + WIN; if ((hi_odd & 1) == 0) hi_odd--;
        n_odd = (u32)((hi_odd - lo_odd) / 2 + 1);
        comp.assign(n_odd, 0);
    }

    void init_first(u64 S, u64 maxN, u64 BS_) {
        BS = BS_;
        set_block(S, maxN);
        size_t np = base_primes.size();
        d.resize(np);
        bsp.resize(np);
        for (size_t i = 1; i < np; i++) {           // skip p=2
            u64 p = base_primes[i];
            u64 base = lo_odd;
            u64 p2 = p * p;
            if (p2 > base) base = p2;               // start from p^2
            u64 m0 = ((base + p - 1) / p) * p;
            if ((m0 & 1) == 0) m0 += p;             // first ODD multiple of p >= base
            d[i] = m0 - lo_odd;                     // even, may be >= 2p
            bsp[i] = (u32)(BS % (2 * p));
        }
    }

    void mark() {
        size_t np = base_primes.size();
        u32 n2 = (u32)(n_odd * 2);
        for (size_t i = 1; i < np; i++) {
            u64 p = base_primes[i];
            u64 dd = d[i];
            if (dd < n2) {
                u32 idx = (u32)(dd >> 1);
                for (; idx < n_odd; idx += (u32)p) comp[idx] = 1;
            }
            if (dd >= BS) {
                d[i] = dd - BS;
            } else {
                u64 p2 = 2 * p;
                if (dd < p2) {
                    u32 b = bsp[i];
                    d[i] = (dd >= b) ? (dd - b) : (dd - b + p2);
                } else {
                    d[i] = (dd + p2 - bsp[i]) % p2;
                }
            }
        }
    }

    bool scan_windows(u64 S, u64 maxN, u64& foundN) {
        pr.clear();
        for (u32 i = 0; i < n_odd; i++)
            if (comp[i] == 0) pr.push_back(lo_odd + 2ULL * i);
        size_t n = pr.size();
        u64 fS = 0;
        {
            size_t lo = 0; while (lo < n && pr[lo] < S) lo++;
            size_t hi = lo; while (hi < n && pr[hi] <= S + WIN) hi++;
            fS = (u64)(hi - lo);
        }
        u64 cur_f = fS, N = S;
        size_t i1 = 0, i2 = 0;
        const u64 INF = ~0ULL;
        while (N <= maxN) {
            u64 T1 = INF;
            while (i1 < n && pr[i1] < N) i1++;
            if (i1 < n) T1 = pr[i1];
            u64 T2 = INF;
            while (i2 < n && (__int128)pr[i2] - (__int128)WIN - 1 < (__int128)N) i2++;
            if (i2 < n) T2 = (u64)((__int128)pr[i2] - (__int128)WIN - 1);
            u64 T = min(T1, T2);
            if (T == INF) T = maxN;
            if (T > maxN) T = maxN;
            if (cur_f == 12) { foundN = N; return true; }
            int delta = 0;
            if (T1 == T) { do { delta -= 1; i1++; } while (i1 < n && pr[i1] == T); }
            if (T2 == T) { do { delta += 1; i2++; } while (i2 < n && (__int128)pr[i2] - (__int128)WIN - 1 == (__int128)T); }
            cur_f += (u64)delta;
            N = T + 1;
        }
        return false;
    }
};

// ---------------- verification ----------------
bool verify(u64 N) {
    int c = 0;
    for (u64 x = N; x <= N + WIN; x++)
        if (is_prime_mr(x)) c++;
    return c == 12;
}

// ---------------- parallel driver ----------------
static atomic<bool> done(false);
static mutex mtx;
static u64 resultN = 0;
static atomic<u64> total_scanned(0);
static atomic<int> active_workers(0);

void worker(u64 start, u64 end, u64 BS) {
    active_workers.fetch_add(1);
    SegSieve seg;
    seg.init_first(start, min(start + BS - 1, end), BS);
    seg.mark();
    u64 S = start;
    while (!done.load()) {
        u64 maxN = min(S + BS - 1, end);
        u64 fN = 0;
        if (seg.scan_windows(S, maxN, fN)) {
            lock_guard<mutex> lk(mtx);
            if (!done.load()) { resultN = fN; done.store(true); }
            active_workers.fetch_sub(1);
            return;
        }
        total_scanned.fetch_add(maxN - S + 1);
        S += BS;
        if (S > end) break;
        seg.set_block(S, min(S + BS - 1, end));
        seg.mark();
    }
    active_workers.fetch_sub(1);
}

int main(int argc, char** argv) {
    // usage: scan_range <start> <end> [bs]
    u64 BS = 20000000ULL;
    if (argc > 1) L = strtoull(argv[1], nullptr, 10);
    if (argc > 2) R = strtoull(argv[2], nullptr, 10);
    if (argc > 3) BS = strtoull(argv[3], nullptr, 10);
    if (L > R) { cerr << "invalid range" << endl; return 1; }
    cerr << "range [" << L << ", " << R << "] bs=" << BS << endl;
    cerr << "Generating base primes up to " << SMAX << " ..." << endl;
    auto t0 = chrono::steady_clock::now();
    gen_base_primes();
    auto t1 = chrono::steady_clock::now();
    cerr << "  base primes: " << base_primes.size() << " in "
         << chrono::duration<double>(t1 - t0).count() << " s" << endl;

    cerr << "Scanning [" << L << ", " << R << "] block=" << BS << " ..." << endl;

    u64 span = R - L + 1;
    u64 mid = L + span / 2;
    auto tp = chrono::steady_clock::now();
    thread t1w(worker, L, mid, BS);
    thread t2w(worker, mid + 1, R, BS);

    while (!done.load() && active_workers.load() > 0) {
        this_thread::sleep_for(chrono::seconds(20));
        if (done.load()) break;
        double el = chrono::duration<double>(chrono::steady_clock::now() - tp).count();
        u64 sc = total_scanned.load();
        double rate = (double)sc / el / 1e6;
        cerr << "  elapsed " << (int)el << " s, scanned " << sc
             << " (" << rate << " Mpos/s, " << (double)sc / (double)(R - L + 1) * 100.0 << "% of range)" << endl;
    }
    t1w.join();
    t2w.join();
    auto t2 = chrono::steady_clock::now();
    cerr << "Scan finished in " << chrono::duration<double>(t2 - tp).count() << " s" << endl;

    if (resultN) {
        bool ok = verify(resultN);
        cout << "FOUND N = " << resultN << endl;
        cout << "  exactly " << WIN + 1 << " numbers [N, N+" << WIN << "] have 12 primes: "
             << (ok ? "TRUE" : "FALSE") << endl;
        cout << "  N in [" << L << ", " << R << "]: " << (resultN >= L && resultN <= R ? "TRUE" : "FALSE") << endl;
        if (ok && resultN >= L && resultN <= R) {
            vector<u64> ps;
            for (u64 x = resultN; x <= resultN + WIN; x++)
                if (is_prime_mr(x)) ps.push_back(x);
            cout << "  primes: ";
            for (size_t i = 0; i < ps.size(); i++) cout << (i ? ", " : "") << ps[i];
            cout << endl;
        }
    } else {
        cout << "NO SOLUTION FOUND in range" << endl;
    }
    return 0;
}
