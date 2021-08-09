
function water_rel_perm(sw::T, swr::P, sor::P, krw0::P, nw::P) where {T,P <: Float64}

    S  = (sw - swr) / (1 - sor - swr) 

    if swr <= sw <= 1 - sor
        krw = krw0 * S^nw
    elseif 1 - sor < sw <= 1.0
        krw = krw0
    elseif sw < swr
        krw = 0.0
    end
    
end


function water_rel_perm(sw::AbstractVector{T}, swr::P, sor::P, krw0::P, nw::P) where {T,P <: Float64}
    water_rel_perm.(sw, swr, sor, krw0, nw)
end


function oil_rel_perm(sw::T, swr::P, sor::P, kro0::P, no::P) where {T,P <: Float64}

    S  = (sw - swr) / (1 - sor - swr) 

    if swr <= sw <= 1 - sor
        kro = kro0 * (1 - S)^no
    elseif 1 - sor < sw <= 1.0
        kro = 0.0
    elseif sw < swr
        kro = kro0
    end
end


function oil_rel_perm(sw::AbstractVector{T}, swr::P, sor::P, kro0::P, no::P) where {T,P <: Float64}
    oil_rel_perm.(sw, swr, sor, kro0, no)
end


function oil_rel_perm(sw::Union{T, AbstractArray{T}}, kr::RelPerms) where {T <: Float64}
    swr = kr.swr
    sor = kr.sor
    kro0 = kr.kro0
    no = kr.no
    
    oil_rel_perm.(sw, swr, sor, kro0, no)
end


function water_rel_perm(sw::Union{T, AbstractArray{T}}, kr::RelPerms) where {T <: Float64}
    swr = kr.swr
    sor = kr.sor
    krw0 = kr.krw0
    nw = kr.nw
    
    water_rel_perm.(sw, swr, sor, krw0, nw)
end



function kro_derivative(sw::T, swr::P, sor::P, kro0::P, no::P) where {T,P <: Float64}
    derivative.(sw ->  oil_rel_perm(sw, swr, sor, kro0, no), sw)
end


function krw_derivative(sw::AbstractVector{T}, swr::P, sor::P, krw0::P, nw::P) where {T,P <: Float64}
    derivative.(sw ->  water_rel_perm(sw, swr, sor, krw0, nw), sw)
end

function fractional_flow(sw::T, swr::P, sor::P, krw0::P, kro0::P, nw::P, no::P, μw::P, μo::P) where {T,P <: Float64}
    krw = water_rel_perm(sw, swr, sor, krw0, nw)
    kro = oil_rel_perm(sw, swr, sor, kro0, no)
    λw = krw / μw
    λo = kro / μo

    λw ./ (λw .+ λo)
end


function fractional_flow(sw::AbstractVector{T}, swr::P, sor::P, krw0::P, kro0::P, nw::P, no::P, μw::P, μo::P) where {T,P <: Float64}
    fractional_flow.(sw, swr, sor, krw0, kro0, nw, no, μw, μo)
end


function fractional_flow(sw::Union{T,AbstractVector{T}}, kr::RelPerms, μw::P, μo::P) where {T,P <: Float64}
    swr = kr.swr
    sor = kr.sor
    krw0 = kr.krw0
    nw = kr.nw
    kro0 = kr.kro0
    no = kr.no
    
    krw = water_rel_perm(sw, swr, sor, krw0, nw)
    kro = oil_rel_perm(sw, swr, sor, kro0, no)
    λw = krw / μw
    λo = kro / μo

    λw ./ (λw .+ λo)
end



function fw_derivative(sw::Union{T,AbstractVector{T}}, swr::P, sor::P, krw0::P, kro0::P, nw::P, no::P, μw::P, μo::P) where {T,P <: Float64}
    derivative.(sw ->  fractional_flow(sw, swr, sor, krw0, kro0, nw, no, μw, μo), sw)
end

function fw_derivative(sw::Union{T,AbstractVector{T}}, kr::RelPerms, μw::P, μo::P) where {T,P <: Float64}
    swr = kr.swr
    sor = kr.sor
    krw0 = kr.krw0
    nw = kr.nw
    kro0 = kr.kro0
    no = kr.no
    
    derivative.(sw ->  fractional_flow(sw, swr, sor, krw0, kro0, nw, no, μw, μo), sw)
end
