function plotPath2d(roi::Vector{Vector{Float64}})
    plt.figure()
    plt.scatter(roi[1][1], roi[1][2])
    show()
end