module CharacterHistoriesMakieExt

import Makie, CharacterHistories

function CharacterHistories.treeplot(
    tree::CharacterHistories.Root,
    model::CharacterHistories.Mk
)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1], xreversed = true)

    Makie.hidespines!(ax)
    Makie.hidedecorations!(ax)

    x, y, states = CharacterHistories.coordinates(tree)

    c_states = String[]
    points = Makie.Point2{Float64}[]
    for index in keys(x)
        for i in 2:length(x[index])
            from = Makie.Point(x[index][i-1], y[index])
            to   = Makie.Point(x[index][i],   y[index])
        
            push!(points, from)
            push!(points, to)
            push!(c_states, states[index][i-1])
        end
    end
    
    #state_space = unique(c_states)
    state_space = model.state_space
    tbl = Dict(x => i for (i,x) in enumerate(state_space))
    cs = [tbl[x] for x in c_states]

    cmap = Makie.cgrad([:black, :red, :green], categorical=true)

    Makie.linesegments!(ax, points, color = cs, colormap = cmap, linewidth = 2)

    return(fig)

end


end