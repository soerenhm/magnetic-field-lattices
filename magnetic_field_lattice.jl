include("core.jl")
import PyPlot as plt; plt.pygui(true)


λ = 1480e-9
θ = deg2rad(10)
ψ = deg2rad(5)

wfun = PlaneWave(λ, [ψ, π/2], [1, 1])
wsh1 = PlaneWave(λ/2, [θ/2, 0], [1, 1])
wsh2 = PlaneWave(λ/2, [θ/2, π], [1, -1])

x = -5λ:0.1λ:5λ
y = copy(x)
z = [0, λ/2, λ]
Efun = sample(wfun, x, y, z)
Esh1 = sample(wsh1, x, y, z)
Esh2 = sample(wsh2, x, y, z)
J = currentdensity(Efun, Esh1 .+ Esh2)[:,:,3]
Jx, Jy, Jz = map(n -> getindex.(J, n), 1:3)
B = magnetic_field(real.(J), x, y)
Bz = getindex.(B, 3)

plt.figure()
plt.subplot(111, aspect="equal")
plt.pcolormesh(x,y,Bz',shading="auto")

I = argmax(Bz)
z = -5λ:0.1λ:5λ
x0, y0 = x[I[1]], y[I[2]]
x = x0-2λ:0.1λ:x0+2λ
y = y0-2λ:0.1λ:y0+2λ
Efun = sample(wfun, x, y, z)
Esh1 = sample(wsh1, x, y, z)
Esh2 = sample(wsh2, x, y, z)
J = currentdensity(Efun, Esh1 .+ Esh2)
B = [getindex.(magnetic_field(real.(J[:,:,i]), x, y),3) for i in eachindex(z)]
Bz = [B[i][size(B[i]).>>1...] for i in eachindex(B)]