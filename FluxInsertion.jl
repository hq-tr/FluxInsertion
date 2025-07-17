include("/home/trung/_qhe-julia/FQH_state_v2.jl")
using .FQH_states
include("/home/trung/_qhe-julia/Misc.jl")
using .MiscRoutine

using Combinatorics

using ArgMacros

function raised(partition::BitVector, index_set::Vector{Int}, No::Int)
    compare = BitVector(dex2bin(index_set, No)) .⊻ BitVector(dex2bin(index_set.+1, No))
    return partition .⊻ compare
end


function all_valid_raise(partition::BitVector, k::Int, Ne::Int, No::Int)
    all_choices = combinations(bin2dex(partition), k)
    all_results = map(c->raised(partition,c, No), all_choices)
    return all_results[count.(all_results).==Ne]
end

function FQH_quasihole_poly_gen!(start_state::FQH_state_mutable, pos::Vector{T} where T <: Number)
    q = length(pos)
    @assert q≥0
    if length(pos)==0
        return
    else
        FQH_quasihole_poly_gen!(start_state,pos[1:end-1])
        w = pos[end]
        println("Inserting flux at $w")

        @time begin
        if w==0
            for vec in start_state.basis
                insert!(vec,1,0)
            end
        else
            for vec in start_state.basis
                push!(vec,0)
            end
            dim = length(start_state.basis)
            Ne  = count(start_state.basis[1])
            No  = length(start_state.basis[1])
            for k in 0:Ne
                for i in 1:dim
                all_raise = all_valid_raise(start_state.basis[i], k, Ne, No)
                all_coef  = (-w/2)^(Ne-k) * start_state.coef[i]*ones(length(all_raise))
                append!(start_state.basis, all_raise)
                append!(start_state.coef, all_coef) 
                end               
            end
            deleteat!(start_state.basis,1:dim)
            deleteat!(start_state.coef,1:dim)

        end
        end
        return
    end
end

function main()
    #@inlinearguments begin
    #@argumentrequired String fname "-f" "--filename"
    println(" ========================== START ============================")

    print("Working on the disk(1) or sphere(2)? ")
    choice = parse(Int, readline())
    @assert choice==1 || choice==2

    println("File name of the jack polynomial ground state:")
    fname = readline()
    state = readwf(fname; mutable=true)
    #println("ground state = ")
    #display(state)

    if length(state.basis) > 0
        println("How many fluxes to add? ")
        nflux = parse(Int, readline())
        loc = ComplexF64[]
        if choice == 1
            println("Input positions X Y for the inserted flux")
            for i in 1:nflux
                print("    Flux #$i: ")
                locread = readline()
                x,y = map(x->parse(Float64, x), split(locread))
                push!(loc, x + y*im)
            end

        elseif choice == 2
            println("Insert fluxes at positions θ = pπ, ϕ = qπ. Input p q:")
            for i in 1:nflux
                print("    Flux #$i: ")
                locread = readline()
                θ,ϕ = map(x->parse(Float64, x), split(locread))
                
                u = cos(π*(1-θ)/2) * exp(π*im*ϕ/2)     # In the formula θ is taken from the South pole. π-θ to flip the pole.
                v = sin(π*(1-θ)/2) * exp(-π*im*ϕ/2)
                
                push!(loc, 2*u/v)
            end
        end

        @time FQH_quasihole_poly_gen!(state,loc)
        collapse!(state)

        println((length(state.basis), length(state.coef)))

        if choice == 1
            disk_normalize!(state)
            printwf(state;fname="$(fname)_flux_dsk")
        elseif choice==2
            sphere_normalize!(state)
            printwf(state;fname="$(fname)_flux_sph")
        end
    end
end

main()
