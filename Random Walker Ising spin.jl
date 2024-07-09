#Importing packages
using Plots
using StatsBase
using StatsPlots
using Distributions
using LinearAlgebra
using LaTeXStrings
using CSV, DataFrames
using Random


#=Function to impose periodic boundary conditions:
if index site is less than 1, it returns n (the index of the opposite site of the lattice) and vice versa. 
this way, the whole lattice is wrapped around itself, satisfying PBC
=#
function boundary(i,n)
    if i<1
        return n
    elseif i>n
        return 1
    else
        return i
    end
end    


function init_st(n, p)
    L=zeros(n,n); #creating a matrix which will be our lattice

    #assigning spins to lattice site according to probability p 
    for i =1:n
        for j = 1:n
            a=rand()
            if a<p
                L[i,j]=1
            else
                L[i,j]=-1
            end
        end
    end
    L
end


#dimension
n=100
#running for 100 trials
#L= init_st(n, 0.5);

X=zeros(n,n)
MCS=n*n
step=4000
#running each MC step
for temp in [2.27]
    Avg=[]
    for trial =1:1
        L= init_st(n, 0.5);
        X=zeros(n,n)
        for m =1:step*MCS
            b=1/temp
            #choosing a random lattice point
            i=rand(1:n)
            j=rand(1:n)
            S= -2*L[i,j]*(L[boundary(i+1,n),j]+L[boundary(i-1,n),j]+L[i,boundary(j+1,n)]+L[i,boundary(j-1,n)]) #calculating energy difference between the present state and if spin is flipped
            E_diff=-S
            #using the spin flip rules
            if E_diff<=0
                    L[i,j]*=-1
            else
                A=exp(-b*E_diff)
                cnt=rand()
                if cnt<A
                    L[i,j]*=-1
                end
            end
            if m%MCS==0
                X=X+L
            end
            if m in [step*MCS]
                #density!(reshape(X,(n*n,1)), label="T=$para")
                append!(Avg,reshape(X,(n*n,1)))
            end
        end
        df = DataFrame(trials=vec(reshape(X,(n*n,1))))
        CSV.write("Random_Walk distributionT2.27_$step MCS.csv", df) 
    end  
end  


#dimension
n=100
#plot()
#Metropolis
#defining range of beta for which to run
para=0
#running for 100 trials
#L= init_st(n, 0.5);

X=zeros(n,n)
MCS=n*n
#running each MC step
#change temperature value accordingly
for temp in [2.27]
    Avg=[]
    for trial =1:500
        L= init_st(n, 0.5);
        X=zeros(n,n)
        for m =1:step*MCS
            b=1/temp
            #choosing a random lattice point
            i=rand(1:n)
            j=rand(1:n)
            S= -2*L[i,j]*(L[boundary(i+1,n),j]+L[boundary(i-1,n),j]+L[i,boundary(j+1,n)]+L[i,boundary(j-1,n)]) #calculating energy difference between the present state and if spin is flipped
            E_diff=-S
            #using the spin flip rules
            if E_diff<=0
                    L[i,j]*=-1
            else
                A=exp(-b*E_diff)
                cnt=rand()
                if cnt<A
                    L[i,j]*=-1
                end
            end
            if m%MCS==0
                X=X+L
            end
            if m in [step*MCS]
                #density!(reshape(X,(n*n,1)), label="T=$para")
                append!(Avg,reshape(X,(n*n,1)))
            end
        end
        col= CSV.read("Random_Walk distributionT2.27_$step MCS.csv", DataFrame)
        df = DataFrame(trials=vec(reshape(X,(n*n,1))))
        col1=hcat(col, df, makeunique=true)
        CSV.write("Random_Walk distributionT2.27_$step MCS.csv", col1) 
    end
end

        
	

