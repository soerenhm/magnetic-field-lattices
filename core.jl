using StaticArrays


function rot(θ, ϕ)
  cosθ, sinθ = cos(θ), sin(θ)
  cosϕ, sinϕ = cos(ϕ), sin(ϕ)
  Ry = @SMatrix[cosθ 0 sinθ; 0 1 0; -sinθ 0 cosθ]
  Rz = @SMatrix[cosϕ -sinϕ 0; sinϕ cosϕ 0; 0 0 1]
  Rz * Ry
end
rotdeg(θ, ϕ) = rot(deg2rad(θ), deg2rad(ϕ))


struct JonesVector
  x::ComplexF64
  y::ComplexF64
  function JonesVector(x::ComplexF64, y::ComplexF64)
    r = sqrt(x*x' + y*y')
    iszero(r) && throw(ArgumentError("Vector is zero and not normalizable!"))
    new(x/r, y/r)
  end
end
JonesVector(x::Number, y::Number) = JonesVector(ComplexF64(x), ComplexF64(y))
Base.vec(e::JonesVector) = @SVector[e.x, e.y, 0]


struct UnitVector
  θ::Float64
  ϕ::Float64
  UnitVector(θ::Float64, ϕ::Float64) = new(θ % π, ϕ % 2π)
end
UnitVector(θ::Real, ϕ::Real) = UnitVector(float(θ), float(ϕ))
Base.vec(u::UnitVector) = @SVector[sin(u.θ)*cos(u.ϕ), sin(u.θ)*sin(u.ϕ), cos(u.θ)]
rot(u::UnitVector) = rot(u.θ, u.ϕ)


struct PlaneWave
  λ::Float64
  dir::UnitVector
  pol::JonesVector
  function PlaneWave(λ::Float64, dir::UnitVector, pol::JonesVector)
    λ > 0 || throw(ArgumentError("λ must be positive!"))
    new(λ, dir, pol)
  end
end
PlaneWave(λ::Real, angles::AbstractArray, pol::AbstractArray) = PlaneWave(float(λ), UnitVector(angles...), JonesVector(pol...))
wavenumber(p::PlaneWave) = 2π/p.λ
wavevector(p::PlaneWave) = vec(p.dir) * wavenumber(p)
polarization(p::PlaneWave) = rot(p.dir) * vec(p.pol)

evaluate_planewave(r, k, e) = e .* exp(1im*k'*r)
evaluate(pw::PlaneWave, r) = evaluate_planewave(r, wavevector(pw), polarization(pw))

function sample(pw::PlaneWave, x, y, z)
  E = zeros(SVector{3,ComplexF64}, length.((x,y,z)))
  sample!(E, wavevector(pw), polarization(pw), x, y, z)
end
function sample!(E, k, e, x, y, z)
  for I in CartesianIndices(E)
    a, b, c = Tuple(I)
    r = @SVector[x[a],y[b],z[c]]
    E[I] = evaluate_planewave(r, k, e)
  end
  E
end


# Assumes isotropic η...
function currentdensity(Eω, E2ω)
  J = similar(Eω)
  currentdensity!(J, Eω, E2ω)
end
function currentdensity!(J, Eω, E2ω)
  for I in CartesianIndices(J)
    e1, e2 = Eω[I], E2ω[I]
    J[I] = (e1 .* e1) .* conj(e2)
  end
  J
end


function magnetic_field(J, x, y)
  B = similar(J)
  magnetic_field!(B, J, x, y)
end
function magnetic_field!(B, J, x, y)
  dA = step(x) * step(y)
  for I in CartesianIndices(B)
    B[I] = _magnetic_field(J, I[1], I[2], x, y)*dA*1e-7
  end
  B
end

function magnetic_field(J, x, y, z)
  B = similar(J)
  magnetic_field!(B, J, x, y, z)
end
function magnetic_field!(B, J, x, y, z)
  for I in CartesianIndices(B)
    B[I] = _magnetic_field(J, I, x, y, z)*1e-7
  end
end

function _magnetic_field(J, I, x, y, z)
  T = eltype(J)
  i, j, k = Tuple(I)
  xi, yj, zk = x[i], y[j], z[k]
  Bx, By, Bz = zero(T), zero(T), zero(T)
  for K in CartesianIndices(J)
    K == CartesianIndex(i,j,k) && continue
    Δx, Δy, Δz = x[K[1]]-xi, y[K[2]]-yj, z[K[3]]-zk
    ξ = sqrt(Δx^2 + Δy^2 + Δz^2)
    Bx += (Jy*Δz - Jz*Δy)/ξ
    By += (Jz*Δx - Jx*Δz)/ξ
    Bz += (Jx*Δy - Jy*Δx)/ξ
  end
  @SVector[Bx, By, Bz]
end
function _magnetic_field(J::AbstractMatrix, i, j, x, y)
  T = eltype(eltype(J))
  Bx, By, Bz = zero(T), zero(T), zero(T)
  xi, yj = x[i], y[j]
  for (l, yl) in enumerate(y)
    for (k, xk) in enumerate(x)
      l == j && k == i && continue
      Jx, Jy, Jz = J[k,l]
      Δx = xk - xi
      Δy = yl - yj
      ξ = sqrt(Δx^2 + Δy^2)^3
      Bx += -Jz*Δy/ξ
      By += Jz*Δx/ξ
      Bz += (Jx*Δy - Jy*Δx)/ξ
    end
  end
  @SVector[Bx, By, Bz]
end