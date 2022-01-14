using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using BS_DAE

τ_P = 5.0
τ_Q = 8.0
K_P = 5
K_Q = 0.1
V_r = 1

Aᵤ = (V_r + 2 * K_Q ) / (2 * τ_Q * V_r) # * Q
Bᵤ = 1im 
Cᵤ = - 1 / (2 * τ_Q * V_r^2)
Gᵤ = 0
Hᵤ = - K_Q / (τ_Q * V_r)
Mₓ = τ_P
Aₓ = K_P #* P 
Bₓ = - 1 
Cₓ = 0
Gₓ = - K_P 
Hₓ = 0