

export ancestral_character_map


export stochastic_character_map 

function stochastic_character_map(
    tree::Root, 
    model::CharacterHistory
    )
    S = ancestral_state_probabilities(tree, model)

    tree2 = deepcopy(tree)

    redraw_nodes_recursive!(tree2, model, S)

    redraw_branches_recursive!(tree2, model)

    return(tree2)
end

