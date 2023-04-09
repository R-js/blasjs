
// from https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Binary_numeral_system_(base_2)
export default function isqrt(n: number): number {
    // Xₙ₊₁
    let x = n;

    // cₙ
    let c = 0;

    // dₙ which starts at the highest power of four <= n
    let d = 1 << 30; // The second-to-top bit is set.
    // Same as ((unsigned) INT32_MAX + 1) / 2.
    while (d > n)
        d >>= 2;

    // for dₙ … d₀
    while (d != 0) {
        if (x >= c + d) {      // if Xₘ₊₁ ≥ Yₘ then aₘ = 2ᵐ
            x -= c + d;        // Xₘ = Xₘ₊₁ - Yₘ
            c = (c >> 1) + d;  // cₘ₋₁ = cₘ/2 + dₘ (aₘ is 2ᵐ)
        }
        else {
            c >>= 1;           // cₘ₋₁ = cₘ/2      (aₘ is 0)
        }
        d >>= 2;               // dₘ₋₁ = dₘ/4
    }
    return c;                  // c₋₁
}