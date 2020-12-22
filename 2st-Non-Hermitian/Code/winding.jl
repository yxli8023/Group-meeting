# Import package which is necessary
#import Pkg
#Pkg.add("LinearAlgebra")
#Pkg.add("BenchmarkTools")
#Pkg.add("PyPlot")
using LinearAlgebra,PyPlot,DelimitedFiles
# =======================================================
function pauli()
    # 构建Pauli矩阵
    sigx = zeros(Float64,2,2)
    sigy = zeros(ComplexF64,2,2)
    sigz = zeros(Float64,2,2)
    sigm = zeros(ComplexF64,2,2)
    sigp = zeros(ComplexF64,2,2)
    #----
    sigx[1,2] = 1.0
    sigx[2,1] = 1.0
    #----
    sigy[1,2] = -1im
    sigy[2,1] = 1im
    #----
    sigz[1,1] = 1
    sigz[2,2] = -1
    #----
    sigm = (sigx - 1im*sigy)/2.0
    sigp = (sigx + 1im*sigy)/2.0
    return sigp,sigm
end
# =========================================================
function hamset(k::Number,tv::Float64)::Matrix{ComplexF64}
    # 哈密顿量构建
    t1::Float64 = tv
    t2::Float64 = 1.0
    gam::Float64 = 4.0/3.0
    r::Float64 = sqrt(abs((t1 - gam/2)/(t1 + gam/2)))
    sigp,sigm = pauli()
    ham = zeros(ComplexF64,2,2)
    ham = (t1 - gam/2.0 + r*exp(1im*k)*t2)*sigm + (t1 + gam/2.0 + 1/r*exp(-1im*k)*t2)*sigp
    return ham
end
# ===========================================================
function Qmat(tt::Float64,tv::Float64)
    occ::Int64 = 1
    sigz = zeros(Float64,2,2)
    sigz[1,1] = 1
    sigz[2,2] = -1
    h1 = hamset(tt,tv)
    vecR = eigvecs(h1)
    vecL = inv(vecR)'
    q1 = vecR[:,occ]
    q11 = sigz*vecR[:,occ]
    q2 = vecL[:,occ]
    q22 = sigz*vecL[:,occ]
    Q = q11*q22' - q1*q2'
    return Q
end
# ================================================================
function winding(tv::Float64)
    re1::ComplexF64 = 0 + 0im
    kn::Int64 = 150
    qmlist = []
    qlist = []
    dq = []
    for k in 0:kn
        k = 2*pi/kn*k
        Q = Qmat(k,tv)
        append!(qlist,Q[1,2])
        append!(qmlist,Q[2,1])
    end
    #--------------------------
    dq = qlist[2:end] - qlist[1:end-1]
    for k in 1:length(qmlist)-1
        re1 += qmlist[k]*dq[k]
    end
    return re1*1im/(2*pi)
end
# ==================================================================
function main()
    tlist = []
    wlist = []
    for tv in -3:0.01:3
        append!(tlist,tv)
        append!(wlist,winding(tv))
    end
    wlist = map(real,wlist)
    #wlist = map(floor,wlist)
    plot(tlist,wlist)
    xlabel("t1")
    ylabel("Winding")
    savefig("TopoInv.png",bbox_inches="tight")
end
# ===============================
# 非厄密winding number计算
main()