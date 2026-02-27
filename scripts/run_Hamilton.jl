using TDGL_Polycrystal
using HDF5
using CUDA
using Adapt

data = [1, 2, 3, 4, 5]
data_GPU = adapt(CUDABackend(), data)
data_GPU = map(x -> x^2, data_GPU)
println(typeof(data_GPU))
println(data_GPU)
