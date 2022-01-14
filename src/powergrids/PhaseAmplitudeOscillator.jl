using NetworkDynamics: ODEVertex
using PowerDynamics
import PowerDynamics: dimension, symbolsof, construct_vertex 
import PowerDynamics: showdefinition

# Source: https://github.com/PIK-ICoNe/NormalFormPaper
@DynamicNode SchifferApprox(τ_P, τ_Q, K_P, K_Q, V_r, P, Q) begin
    Aᵤ = (V_r + 2 * K_Q * Q) / (2 * τ_Q * V_r)
    Bᵤ = 1im 
    Cᵤ = - 1 / (2 * τ_Q * V_r^2)
    Gᵤ = 0
    Hᵤ = - K_Q / (τ_Q * V_r)
    Mₓ = τ_P
    Aₓ = K_P * P 
    Bₓ = - 1 
    Cₓ = 0
    Gₓ = - K_P 
    Hₓ = 0
    @assert isreal(Mₓ) && Mₓ > 0
    @assert isreal(Aₓ)
    @assert isreal(Bₓ)
    @assert isreal(Cₓ)
    @assert isreal(Gₓ)
    @assert isreal(Hₓ)
end [[ω, dω]] begin
    s = u * conj(i)
    v2 = abs2(u)
    dω = ( Aₓ + Bₓ * ω + Cₓ * v2 + Gₓ * real(s) + Hₓ * imag(s) ) / Mₓ
    du = ( Aᵤ + Bᵤ * ω + Cᵤ * v2 + Gᵤ * real(s) + Hᵤ * imag(s) ) * u
end
