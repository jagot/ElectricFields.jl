compute_rotation(::Nothing) = I
compute_rotation(I::UniformScaling) = I/opnorm(I)

function compute_rotation((ϕ,u)::Tuple{<:Real,<:AbstractVector{<:Real}})
    @assert length(u) == 3
    u = SVector{3}(normalize(u))
    s,c = sincos(ϕ)
    # https://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication
    ucross = SMatrix{3,3}([0    -u[3]  u[2]
                           u[3]  0    -u[1]
                           -u[2] u[1]  0])
    # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    R = c*I + s*ucross + (1-c)*(u*u')
    R/opnorm(R)
end

function compute_rotation(R::AbstractMatrix{T}) where {T<:Real}
    size(R) == (3,3) ||
        throw(DimensionMismatch("Rotation matrix must have dimensions 3×3"))
    rank(R) < 3 &&
        throw(ArgumentError("Rotation matrix singular"))
    R = Matrix{T}(R)
    for i = 2:3
        for j = 1:i-1
            R[:,i] -= dot(R[:,j],R[:,i])*R[:,j]
        end
    end
    SMatrix{3,3}(R/opnorm(R))
end

rotation_angle(R::AbstractMatrix) = acos(clamp((tr(R) - 1)/2, -1, 1))

# https://en.wikipedia.org/wiki/Rotation_matrix#Determining_the_axis
function rotation_axis(R::AbstractMatrix)
    ee = eigen(R)
    i = argmin(abs.(ee.values .- 1))
    normalize(real(ee.vectors[:,i]))
end

export rotation_angle, rotation_axis, rotation_matrix
