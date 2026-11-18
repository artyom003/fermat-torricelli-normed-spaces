using LinearAlgebra, Plots

A = [0.0, 2.0]
B = [-1.0, 0.0]
C = [2.0, 0.0]
pts = [A, B, C]

f(T) = sum(norm(T .- P) for P in pts)

function torricelli_point(pts; tol=1e-10, max_iter=1000)
    T = sum(pts) / length(pts)
    for _ in 1:max_iter
        weights = [1.0 / norm(T .- P) for P in pts]
        T_new = sum(w .* P for (w, P) in zip(weights, pts)) / sum(weights)
        if norm(T_new - T) < tol
            return T_new
        end
        T = T_new
    end
    return T
end

T = torricelli_point(pts)
centroid = sum(pts) / 3

println("Տորիչելիի կետը T = ", round.(T, digits=4))
println("Հեռավորությունների գումարը T-ում = ", round(f(T), digits=4))
println("Ծանրության կենտրոնում գումարը     = ", round(f(centroid), digits=4))

angle_deg(u, v) = acosd(clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0))
println("∠ATB = ", round(angle_deg(A - T, B - T), digits=2), "°")
println("∠BTC = ", round(angle_deg(B - T, C - T), digits=2), "°")
println("∠CTA = ", round(angle_deg(C - T, A - T), digits=2), "°")

p = plot([A[1], B[1], C[1], A[1]], [A[2], B[2], C[2], A[2]],
         lw=2, color=:black, label="△ABC", aspect_ratio=:equal,
         title="Դասական Ֆերմա-Տորիչելիի խնդիրը", xlabel="x", ylabel="y")

for P in pts
    plot!(p, [T[1], P[1]], [T[2], P[2]], lw=2, ls=:dash, color=:blue, label="")
end

scatter!(p, [A[1], B[1], C[1]], [A[2], B[2], C[2]], ms=7, color=:red, label="A, B, C")
scatter!(p, [T[1]], [T[2]], ms=8, color=:green, label="Տորիչելիի կետ (T)")
scatter!(p, [centroid[1]], [centroid[2]], ms=6, color=:orange, shape=:diamond, label="Ծանրության կենտրոն")

savefig(p, "classical_fermat_torricelli.png")
display(p)