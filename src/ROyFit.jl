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


"""
    function write_prefit_summary(fname::AbstractString, detx_fname::AbstractString)

Perform the prefit on a given file and create a CSV summary file.
"""
function write_prefit_summary(fname::AbstractString, detx_fname::AbstractString)
    run_id = mc_run_id(fname)
    calib = read_calibration(detx_fname)
    reader = KM3NeT.EventReader(fname, load_tracks=true)
    output_fname = basename(fname) * ".royprefit.csv"
    out_fobj = open(output_fname, "w")
    write(out_fobj, "run_id,event_id,n_muons,muon_energy,dir_x,dir_y,dir_z,alpha\n")
    for (event_id, event) in enumerate(reader)
        n_muons = length([t for t in event.mc_tracks if t.particle_type == 5])
        muon_energy = sum([t.E for t in event.mc_tracks if t.particle_type == 5])
        write(out_fobj, "$run_id,$event_id,$n_muons,$muon_energy,")
        reco = ROyFit.prefit(event.hits, calib)
        if isa(reco, NoRecoTrack)
            write(out_fobj, "nan,nan,nan,nan\n")
            continue
        end
        muon = event.mc_tracks[2]
        alpha = angle(Direction(reco.dir), muon.dir)
        write(out_fobj, "$(reco.dir.x),$(reco.dir.y),t$(reco.dir.z),$alpha\n")
    end
    close(out_fobj)
end


end # module
