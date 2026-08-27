include("/home/trung/_qhe-julia/FQH_state_v2.jl")
using .FQH_states
include("/home/trung/_qhe-julia/Misc.jl")
using .MiscRoutine

using Combinatorics
using LinearAlgebra

using ArgMacros

const tol = 1e-14
function raised(partition::BitVector, index_set::Vector{Int}, No::Int)
    compare = BitVector(dex2bin(index_set, No)) .⊻ BitVector(dex2bin(index_set.+1, No))
    return partition .⊻ compare
end


function all_valid_raise(partition::BitVector, k::Int, Ne::Int, No::Int)
    all_choices = combinations(bin2dex(partition), k)
    all_results = map(c->raised(partition,c, No), all_choices)
    return all_results[count.(all_results).==Ne]
end

function FQH_quasihole_poly_gen!(start_state::FQH_state_mutable, pos::Vector{T} where T <: Number;quiet=false)
    q = length(pos)
    @assert q≥0
    if length(pos)==0
        return
    else
        FQH_quasihole_poly_gen!(start_state,pos[1:end-1])
        w = pos[end]
        if !quiet
            println("Inserting flux at $w")
        end

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
        end # end of @time
        return
    end
end

function main()

    @inlinearguments begin
        @argumentrequired String fname "-f" "--filename"
        @argumentrequired String geom "-g" "--geometry"
        @argumentrequired Int nflux "-n" "--num-flux"
        @argumentflag overwrite "--overwrite"
        @argumentflag saveraw "--save-raw"
    end

    if !isfile(fname)
        println("File '$(fname)' not found. Terminated")
        return
    end

    geom = lowercase(geom)
    @assert(geom in ["sphere","disk"], "'geometry' argument must be either 'sphere' or 'disk'.")

    outputname = "$(fname)_add_$(nflux)_flux"
    if isdir(outputname)
        if length(readdir(outputname)) > 0
            if !overwrite
                println("Non-empty output directory '$(outputname)' exists.")
                println("Remove the directory, first, or run the program with flag '--overwrite' to overwrite directory.")
                return
            else
                rm(outputname,recursive=true)
                mkdir(outputname)
            end
        end
    else
        mkdir(outputname)
    end


    ground_state = readwf(fname)
    Ne = count(ground_state.basis[1])
    println("$Ne electrons.")

    guess_num_state = binomial(Ne + nflux,Ne)
    guess_num_state += guess_num_state < 5 ? guess_num_state : 5 # Add an overcounting


    vortex_states = FQH_state_mutable[]
    for i in 1:guess_num_state
        print("\rGenerating state $i out of $guess_num_state      ")
        #state = FQH_state_mutable(copy(ground_state.basis),copy(ground_state.coef)) # create a mutable copy of the ground state
        state = readwf(fname;mutable=true)
        loc = rand(ComplexF64,nflux) # nflux random complex numbers
        FQH_quasihole_poly_gen!(state,loc;quiet=true)
        collapse!(state)
        if geom == "disk"
            disk_normalize!(state)
        elseif geom == "sphere"
            sphere_normalize!(state)
        end
        push!(vortex_states,state)

        if saveraw
            printwf(state;fname="$outputname/raw_$(i-1)")
        end
        state = nothing
        GC.gc()
    end
    println("Done")

    # Collate and orthonormalize
    all_basis, all_coef = collate_many_vectors(vortex_states; separate_out=true, collumn_vector=true)

    #println("Orthonormalizing basis using QR decomposition")
    #@time all_coef_ortho = transpose(Matrix(qr(all_coef).Q))

    println("Orthonormalizing basis using overlap matrix")
    ov_matrix = transpose(conj.(all_coef)) * all_coef
    println(size(ov_matrix))

    E = eigvals(Hermitian(ov_matrix))
    v = eigvecs(Hermitian(ov_matrix))

    all_coef_ortho = [v[:,i] for i in 1:length(E) if E[i] > tol]
    d   = length(all_coef_ortho)

    println("$d orthonormal states")


    for i in 1:d
        print("\rSaving state $i out of $d      ")
        coef = all_coef * all_coef_ortho[i] 
        s = wfnormalize(FQH_state(all_basis,coef))
        println("Check norm = $(wfnorm(s))")
        printwf(s;fname="$(outputname)/vec_$(i-1)")
    end
    println("Done!")


   
end

@time main()
