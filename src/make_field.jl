test_symbol_walk(node, params) = node
test_symbol_walk(node::Symbol, params) = Expr(:call, :in, Expr(:quote, node),
                                              Expr(:call, :keys, params))
test_symbol_walk(node::Expr, params) =
    Expr(node.head, test_symbol_walk.(node.args, params)...)

walk(node, params) = node
walk(node::Symbol, params) = node

get_reference(r::Symbol,params) = Expr(:ref, params, Expr(:quote, r))
get_reference(r,params) = r

get_symbol(s::Symbol, params) = isdefined(s) ? s : get_reference(s, params)
get_symbol(s, params) = s

function walk(node::Expr, params)
    if node.head ∈ [:line, :quote]
        node
    elseif node.head == :if
        args = [test_symbol_walk(node.args[1], params)]
        append!(args, [walk(a, params) for a in node.args[2:end]])
        Expr(node.head, args...)
    else
        args = [walk(a, params) for a in node.args]
        if node.head == :(=)
            r = get_reference(args[1], params)
            if typeof(args[1]) == Symbol
                conv = Expr(:call, :|>, args[2],
                            Expr(:ref, base_units,
                                 Expr(:quote, args[1])))
                Expr(node.head, r, conv)
            else
                Expr(node.head, args...)
            end
        else
            args = get_symbol.(args, params)
            Expr(node.head, args...)
        end
    end
end

macro namespace!(exprs, params)
    local tree = walk(exprs, esc(params))
    quote
        $tree()
    end
end

function test_field_parameters(field_params, set)
    info = set ∩ keys(field_params)
    set_string = join(set, ", ", " and ")

    length(info) == 0 &&
        error("Need to provide one of $(set_string)")
    length(info) > 1 &&
        error("Can only specify one of $(set_string)")

    info
end

function make_field(field_params::Dict{Symbol,Union{Number,Quantity}})
    carrier_info = test_field_parameters(field_params, [:λ, :T, :f, :ν, :ω])
    amplitude_info = test_field_parameters(field_params, [:I₀, :E₀])

    for k in keys(field_params)
        field_params[k] = get_unitful_quantity(field_params, k)
    end

    @namespace!(field_params) do
        if λ || T
            if λ
                T = λ/u"c"
            elseif T
                λ = T/u"c"
            end
            ν = 1/λ
            f = 1/T
            ω = 2π*u"rad"*f
        else # ∝ Frequency specified
            if f
                ν = f/u"c"
                ω = 2π*u"rad"*f
            elseif ν
                f = ν*u"c"
                ω = f/(2π*u"rad")
            elseif ω
                f = ω/(2π*u"rad")
                ν = f/u"c"
            end
            T = 1/f
            λ = 1/ν
        end

        if I₀
            E₀ = √(2I₀/u"ε0*c")
        elseif E₀
            I₀ = u"ε0*c"/2*E₀^2
        end
    end

    field_params
end
