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
