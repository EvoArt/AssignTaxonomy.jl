# R's default random number generator (Mersenne-Twister with R's seeding),
# so that `seed = s` reproduces the bootstrap draws R makes after `set.seed(s)`.
# Port of src/main/RNG.c (R >= 3.6 defaults; runif uses unif_rand directly).

mutable struct RMersenneTwister
    mt::Vector{UInt32}   # 624 state words
    mti::Int
end

"""
    RMersenneTwister(seed)

Generator in the state R is in after `set.seed(seed)` (default kind).
"""
function RMersenneTwister(seed::Integer)
    s = UInt32(mod(seed, 2^32))
    for _ in 1:50                     # initial scrambling
        s = 0x00010dcd * s + 0x00000001
    end
    words = Vector{UInt32}(undef, 625)
    for j in 1:625
        s = 0x00010dcd * s + 0x00000001
        words[j] = s
    end
    # words[1] is R's dummy[0] (mti), reset to N by FixupSeeds
    return RMersenneTwister(words[2:625], 624)
end

function _mt_next!(r::RMersenneTwister)
    N = 624; M = 397
    mt = r.mt
    if r.mti >= N
        @inbounds for kk in 1:N
            y = (mt[kk] & 0x80000000) | (mt[mod1(kk + 1, N)] & 0x7fffffff)
            mt[kk] = mt[mod1(kk + M, N)] ⊻ (y >> 1) ⊻ ((y & 0x1) == 0x1 ? 0x9908b0df : 0x00000000)
        end
        r.mti = 0
    end
    r.mti += 1
    y = mt[r.mti]
    y ⊻= y >> 11
    y ⊻= (y << 7) & 0x9d2c5680
    y ⊻= (y << 15) & 0xefc60000
    y ⊻= y >> 18
    return y
end

const _I2_32M1 = 2.328306437080797e-10

"R's `unif_rand()`: a uniform on the open interval (0, 1)."
function unif_rand(r::RMersenneTwister)
    v = Float64(_mt_next!(r)) * 2.3283064365386963e-10
    v <= 0.0 && return 0.5 * _I2_32M1
    1.0 - v <= 0.0 && return 1.0 - 0.5 * _I2_32M1
    return v
end

"`n` draws as R's `runif(n)`."
runif(r::RMersenneTwister, n::Integer) = [unif_rand(r) for _ in 1:n]
