@testset "Linear Extrapolation" begin
    #test basic functionality
    xdata = [1,2]
    ydata = [2,4]
    m = TDGL_Polycrystal.lin_ext(xdata,ydata)
    @test m ≈ 2

    #test interpolation on log-log
    Jdata = [1e-2,1e-4]
    Edata = [1e-3,1e-5]
    Etarget = 1e-4
    m = TDGL_Polycrystal.lin_ext(log.(Jdata),log.(Edata))
    logJ = TDGL_Polycrystal.invert_linear(m,log(Jdata[end]),log(Edata[end]),log(Etarget))
    J = exp(logJ)
    @test (J < Jdata[1]) & (J > Jdata[end])
end

@testset "Period Average" begin
    #test basic functionality
    data = [1,2,3,4,5]
    period = 2
    avg = TDGL_Polycrystal.period_avg(data, period)
    @test avg == 4

    #test with uneven length data
    data = [1,2,3,4]
    period = 3
    avg = TDGL_Polycrystal.period_avg(data, period)
    @test avg == 2.5
end

@testset "Quaternion" begin
    v = [1.0,0.0,0.0]
    u = [0.0,0.0,1.0]
    q = TDGL_Polycrystal.attitude_quaternion(π/2,v)
    r = TDGL_Polycrystal.rotate(u,q)
    @test r ≈ [0.0,-1.0,0.0] atol=1e-6
end

@testset "Anti-Aliasing" begin
    @testset "Sample" begin
        grid = reshape(1.0:144.0,(12,12))

        avg1 = TDGL_Polycrystal.sample(grid,[2,2],[1,1])
        @test avg1 == 14

        avg2 = TDGL_Polycrystal.sample(grid,[12,12],[1,1])
        @test avg2 == 137.5
    end
    
    @testset "Low Resolution Grid" begin
        grid = reshape(1.0:15.0^2,(15,15))
        newgrid = TDGL_Polycrystal.lower_resolution(grid,3)
        @test size(newgrid) == (5,5)
        @test sum(newgrid)*9 == sum(grid)

        newgrid = TDGL_Polycrystal.lower_resolution(grid,5)
        @test size(newgrid) == (3,3)
        @test sum(newgrid)*25 == sum(grid)
    end
end