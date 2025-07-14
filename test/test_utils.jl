@testset "Linear Extrapolation" begin
    #test basic functionality
    xdata = [1,2]
    ydata = [2,4]
    m = TDGL2D.lin_ext(xdata,ydata)
    @test m ≈ 2

    #test interpolation on log-log
    Jdata = [1e-2,1e-4]
    Edata = [1e-3,1e-5]
    Etarget = 1e-4
    m = TDGL2D.lin_ext(log.(Jdata),log.(Edata))
    logJ = TDGL2D.invert_linear(m,log(Jdata[end]),log(Edata[end]),log(Etarget))
    J = exp(logJ)
    @test (J < Jdata[1]) & (J > Jdata[end])
end

@testset "Quaternions" begin
    #test quaternion creation
    θ = π/4
    q = attitude_quaternion(θ,[1,0,0])
    @test q ≈ [0.9238795325112867, 0.3826834323650898, 0.0, 0.0]

    #test quaternion multiplication
    q1 = attitude_quaternion(θ,[1,0,0])
    q2 = attitude_quaternion(θ,[0,1,0])
    q3 = q1 * q2
    @test q3 ≈ [0.7071067811865476, 0.7071067811865475, 0.0, 0.0]

    #test quaternion rotation
    v = [1.0, 0.0, 0.0]
    rotated_v = rotate(v, q)
    @test rotated_v ≈ [0.7071067811865476, 0.7071067811865475, 0.0] 
end