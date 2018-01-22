module ROyFit

using KM3NeT


const n = 1.3797
const c = 299792458.0 / 1e9


function prefit(hits::Vector{Hit}, calib::Calibration; min_n_doms = 3, tmax = 15, n_hits = 4)
    shits = nfoldhits(hits, tmax, n_hits)
    n_doms = length(Set(sh.dom_id for sh in shits))
    if n_doms < min_n_doms
        return NoRecoTrack()
    end
    chits = calibrate(shits, calib)
    sort!(chits)
    t0 = chits[div(length(shits), 2)].t
    pos, dir = svdfit(matrix([h.pos for h in chits]))
    if (last(chits).pos - first(chits).pos) ⋅ dir  < 0.0
        dir *= -1
    end
    return RecoTrack(dir, pos, t0)
end


struct SingleDUParameters
    t::Vector{HitTime}
    z::Vector{Float32}
    tot::Vector{ToT}
    σ_t::Float32
    mean_tot::Float32

    function SingleDUParameters(t, z, tot, σ_t)
        mean_tot = mean(tot)
        return new(t, z, tot, σ_t, mean_tot)
    end
end


function γ_d(p::SingleDUParameters, uz, zc, dc)
    return n/sqrt(n^2-1) * sqrt(dc^2 + ((p.z-zc)^2) * (1-uz^2))
end


function γ_t(p::SingleDUParameters, uz, zc, dc, tc)
    return tc + ((p.z-zc)*uz+(n^2 - 1) * γ_d(p, uz, zc, dc) / n ) / c
end


function cosθ(p::SingleDUParameters, uz, zc, dc)
    return (1 - uz^2) * (p.z - zc) / γ_D(p, uz, zc, dc) + uz/n
end


function make_quality_function(p::SingleDUParameters)
    function quality_function(uz, zc, dc, tc)
        tot_weights = 2p.tot / cosθ(p, uz, zc, dc) + 1.0
        mean_tot = mean(tot_weights)
        d = sqrt(p.d1^2 + γ_d(p, uz, zc, dc)^2)
        #return sum(γ_t(p, uz, zc, dc, tc) - p.t)^2 / p.σ_t^2 +
    end
    return quality_function
end


end # module
