# Import package which is necessary
#import Pkg
#Pkg.add("LinearAlgebra")
#Pkg.add("BenchmarkTools")
#Pkg.add("PyPlot")
using LinearAlgebra,PyPlot,DelimitedFiles
# ==============================================================
#--------------------------------
function matrixSet(xn::Int64,t1::Float64,t2::Float64,t3::Float64,gam::Float64)
    ham = zeros(ComplexF64,xn*2,xn*2)
    #----------------------------
    sigx = zeros(Float64,2,2)
    sigx[1,2] = 1.0
    sigx[2,1] = 1.0
    sigy = zeros(ComplexF64,2,2)
    sigy[1,2] = -1im
    sigy[2,1] = 1im
    #-----------------------------
    for k in 0:xn-1
        if k == 0    # First line     
           for m1 in 1:2
                for m2 in 1:2
                    ham[m1,m2] = t1*sigx[m1,m2] + 1im*gam/2.0*sigy[m1,m2]
                    ham[m1,m2 + 2] = (t2 + t3)/2.0*sigx[m1,m2] - 1im*(t2 - t3)/2.0*sigy[m1,m2]  
                end
            end
        elseif k == xn-1
            for m1 in 1:2
                for m2 in 1:2
                    ham[k*2 + m1,k*2 + m2] = t1*sigx[m1,m2] + 1im*gam/2.0*sigy[m1,m2]
                    ham[k*2 + m1,k*2 + m2 - 2] = (t2 + t3)/2.0*sigx[m1,m2] + 1im*(t2 - t3)/2.0*sigy[m1,m2]  
                end
            end
        else
            for m1 in 1:2
                for m2 in 1:2
                    ham[k*2 + m1,k*2 + m2] = t1*sigx[m1,m2] + 1im*gam/2.0*sigy[m1,m2]
                    #   right hopping
                    ham[k*2 + m1,k*2 + m2 + 2] = (t2 + t3)/2.0*sigx[m1,m2] - 1im*(t2 - t3)/2.0*sigy[m1,m2]
                    #   left hopping
                    ham[k*2 + m1,k*2 + m2 - 2] = (t2 + t3)/2.0*sigx[m1,m2] + 1im*(t2 - t3)/2.0*sigy[m1,m2]  
                end
            end
        end
    end
    #----------------------------
    return ham
end
# ======================================================================
function eigHam(k::Float64)
    t1 = 1.0 + 0im
    t2 = 1.0 + 0im
    t3 = 0
    gam = 1.0
   #-------------------
    dx = t1 + (t2 + t3)*cos(k)
    dy = (t2 - t3)*sin(k)
    return sqrt(dx^2 + (dy + 1im*gam/2.0)^2)
end
# ========================================================================
function edgeBand(xn::Int64)
    en::Int64 = 100;
    valSet = zeros(ComplexF64,2*en + 1,xn*2)
    valSetre = zeros(ComplexF64,2*en + 1,xn*2)
    valSetim = zeros(ComplexF64,2*en + 1,xn*2)
    tlist = []
    for m1 in -en:en
        t1 = -3.0*m1/en
        append!(tlist,t1)
        ham = matrixSet(xn,t1,1.0,0.0,4.0/3.0)
        val = eigvals(ham)
        repart = map(real,val)
        impart = map(imag,val)
        val = map(abs,val)
        sort!(val)
        sort!(repart)
        sort!(impart)
        valSetre[m1 + en + 1,:] = repart
        valSetim[m1 + en + 1,:] = impart
        valSet[m1 + en + 1,:] = val
    end
    PyPlot.figure()
    PyPlot.plot(tlist,valSet)
    xlabel("t1")
    ylabel("|E|")
    savefig("absE.png",bbox_inches="tight")
    PyPlot.figure()
    PyPlot.plot(tlist,valSetre)
    xlabel("t1")
    ylabel("Re(E)")
    savefig("ReE.png",bbox_inches="tight")
    PyPlot.figure()
    PyPlot.plot(tlist,valSetim)
    xlabel("t1")
    ylabel("Im(E)")
    savefig("ImE.png",bbox_inches="tight")
    PyPlot.show()
end
# ======================================================================
function wave(xn::Int64)
    ham = matrixSet(xn,1.0,1.0,0.0,4.0/4.0)
    PyPlot.figure()
    val,vec = eigen(ham)
    PyPlot.plot(1:xn,map(abs,vec[1:xn,:]))
    PyPlot.show()
end
# =======================================================================
function main()
    edgeBand(40)
    wave(40)
end
# ========================================================================
# 能带计算及非厄米趋肤效应
main()