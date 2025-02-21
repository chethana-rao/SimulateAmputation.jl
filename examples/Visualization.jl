fig2 = Figure(size=(1200,800))

ax1 = LScene(fig2[1,1])
# ax2 = LScene(fig1[1,2])

for i in 1:1:2
    hp1 = poly!(ax1,GeometryBasics.Mesh(VV[i],FF[i]), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
end

# hp5 = scatter!(ax1,p_distal; color=:red,markersize=25)
# hp2 = wireframe!(ax1,GeometryBasics.Mesh(Vs,Eb_low), linewidth=5,color=:red)
# hp4 = dirplot(ax1,Vs[indCurve],NB_base.*mean(edgelengths(Fs,Vs))*3; color=:cyan,linewidth=3,scaleval=1.0,style=:from)

# hp7 = mesh!(ax1,GeometryBasics.Mesh(Vn,Fn), color=:white, shading = FastShading, transparency=false)

fig2