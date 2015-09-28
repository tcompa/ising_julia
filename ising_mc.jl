# program:       ising_mc.jl
# author:        tc
# created:       2015-08-20 -- 8:30 CEST
# last-modified: 2015-08-22 -- 18:30 CEST
# description:   implements a local Metropolis sampling for the two-dimensional
#                Ising model


# construct the neighbors table
function init_nbr(L :: Int)
    N = L^2 :: Int
    nbr = Array{Int, 2}
    nbr = zeros(Int, (N, 4))
    for site = 0:N-1
        nbr[site + 1, 1] = div(site, L) * L + mod(site + 1, L) + 1
        nbr[site + 1, 2] = mod(site + L, N) + 1
        nbr[site + 1, 3] = div(site, L) * L + mod(site - 1, L) + 1
        nbr[site + 1, 4] = mod(site - L, N) + 1
    end
    return nbr[:, :]
end


# runs the algorithm for a given set of parameters
function ising_run(L :: Int,
                   beta :: Float64,
                   num_sweeps :: Int,
                   nbr :: Array{Int, 2})
    # initialize configuration S
    S = Array(Int8, L, L)
    fill!(S, 1)
    # declare a set of useful variables
    N = L ^ 2 :: Int
    dE = 0 :: Int
    tot_M = 0 :: Int
    tot_M_sq = 0 :: Int
    acc = 0 :: Int
    # run loop over sweep/steps (1 sweep = N steps)
    for sweep = 1:num_sweeps
        for step = 1:N
            site = rand(1:N)
            dE = 2 * S[site] * sum(S[nbr[site, :]])
            if rand() < exp(- beta * dE)
                S[site] = - S[site]
                acc += 1
            end
        end
        M = abs(sum(S))
        tot_M += M
        tot_M_sq += M^2
    end
    av_M = (tot_M / num_sweeps) :: Float64
    av_M_sq = (tot_M_sq / num_sweeps) :: Float64
    av_m = (av_M / N) :: Float64
    xi = (beta * (av_M_sq - av_M^2) / N) :: Float64
    acc /= (num_sweeps * N)
    return av_m, xi, acc
end


##############################################################################
# Main
##############################################################################

# choose parameters
L = 16 :: Int
num_sweeps = 20000 :: Int
list_beta = [0.15 : 0.05 : 0.65] :: Array{Float64,1}

# construct neighbor table
nbr = init_nbr(L)

# initialize output files
out_data = open(@sprintf("data_L%d.dat", L), "w")
out_logs = open(@sprintf("logs_L%d.dat", L), "w")
write(out_logs, @sprintf("L=%d - num_sweeps=%d\n", L, num_sweeps))
flush(out_logs)

# run (tic/toq are used for timing)
tic()
for beta = list_beta
    println(beta)
    av_m, xi, acc = ising_run(L, beta, num_sweeps, nbr)
    write(out_data, @sprintf("%.4f %.10f %.10f %.4f\n", beta, av_m, xi, acc))
    flush(out_data)
end
tot_time = toq()
write(out_logs, @sprintf("Elapsed time: %f s\n", tot_time))
close(out_data)
close(out_logs)
