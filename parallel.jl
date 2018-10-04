# choose how many workers processes you want
howmanyprocs = 16

# add the correct amount of processes so that there are 4 worker processes
addprocs(howmanyprocs+1-nprocs())
@assert nworkers()==howmanyprocs

# this function is not designed to be fast
@everywhere function getnthprime(n::Int)
    @assert n >=1
    listofprimes = Array{Int,1}()
    push!(listofprimes,2)
    candidate = 2
    while (length(listofprimes) < n)
        candidate+=1
        prime=true
        for k in listofprimes
            if candidate % k == 0
                prime = false
                break
            end
        end
        if prime
            push!(listofprimes,candidate)
        end
    end
    return listofprimes[end]
end

#create a number of 3000x3000 matrices populated with random input
# M = Matrix{Float64}[rand(3000,3000) for i = 1:nummat];
H = collect(200:400:20000)

# warm up these functions 
# get the primes serially
for i in 1:length(H)
    getnthprime(H[i])
end
# get the primes in parallel
pmap(getnthprime,H);

# get the primes serially
tic()
for i in 1:length(H)
    getnthprime(H[i])
end
serial = toc()

# get the primes in parallel
tic()
pmap(getnthprime,H)
parallel = toc()

println("The parallel runs $(serial/parallel) faster than the serial.")