module CharacterHistoriesMakieExt

import Makie, CharacterHistories

function CharacterHistories.treeplot(
    tree::CharacterHistories.Root
)

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1], xreversed = true)

    Makie.hidespines!(ax)
    Makie.hidedecorations!(ax)

    ## horizontal lines
    x, y, states = CharacterHistories.coordinates_horizontal(tree)

    c_states_horizontal = String[]
    points_horizontal = Makie.Point2{Float64}[]
    for index in keys(x)
        for i in 2:length(x[index])
            from = Makie.Point(x[index][i-1], y[index])
            to   = Makie.Point(x[index][i],   y[index])
        
            push!(points_horizontal, from)
            push!(points_horizontal, to)
            push!(c_states_horizontal, states[index][i-1])
        end
    end

    ## vertical lines
    x, y, states = CharacterHistories.coordinates_vertical(tree)

    c_states_vertical = String[]
    points_vertical = Makie.Point2{Float64}[]
    for index in keys(x)
        #for i in 2:length(x[index])
            from = Makie.Point(x[index], y[index][1])
            to   = Makie.Point(x[index], y[index][2])
        
            push!(points_vertical, from)
            push!(points_vertical, to)
            push!(c_states_vertical, states[index])
        #end
    end
    
    points = vcat(points_horizontal, points_vertical)
    c_states = vcat(c_states_horizontal, c_states_vertical)

    state_space = unique(values(CharacterHistories.tipstates(tree)))

    tbl = Dict(x => i for (i,x) in enumerate(state_space))
    cs = [tbl[x] for x in c_states]

    cmap = Makie.cgrad([:black, :red, :green], categorical=true)

    Makie.linesegments!(ax, points, color = cs, colormap = cmap, linewidth = 2)

    return(fig)

end


end ## end module extension