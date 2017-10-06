module ROyFit

using KM3NeT

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

end # module
