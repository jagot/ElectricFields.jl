base_units = Dict(:λ => u"nm",
                  :I => u"W/cm^2",
                  :E => u"V/m",
                  :τ => u"fs",
                  :T => u"s",
                  :f => u"Hz",
                  :ν => u"cm^-1",
                  :ω => u"rad/s")

function get_unitful_quantity(field_params::Dict{Symbol,Union{Number,Quantity}}, sym::Symbol)
    v = field_params[sym]
    typeof(v) <: Quantity ? v : v*base_units[sym]
end

walk(node, params) = node
walk(node::Symbol, params) = node

get_reference(r::Symbol,params) = Expr(:ref, params, Expr(:quote, r))
get_reference(r,params) = r

get_symbol(s::Symbol, params) = isdefined(s) ? s : get_reference(s, params)
get_symbol(s, params) = s

function walk(node::Expr, params)
    if node.head ∈ [:line, :quote]
        node
    else
        args = [walk(a, params) for a in node.args]
        if node.head == :(=)
            args[1] = get_reference(args[1], params)
            dump(node)
        else
            args = get_symbol.(args, params)
        end
        Expr(node.head, args...)
    end
end

macro namespace!(exprs, params)
    local tree = walk(exprs, params)
    println(tree)
    quote
        $tree
    end
end

function test_field_parameters(field_param_keys, set)
    info = set ∩ field_param_keys
    set_string = join(set, ", ", " and ")

    length(info) == 0 &&
        error("Need to provide one of $(set_string)")
    length(info) > 1 &&
        error("Can only specify one of $(set_string)")

    info
end

function make_field(field_params::Dict{Symbol,Union{Number,Quantity}})
    ks = keys(field_params)
    carrier_info = test_field_parameters(ks, [:λ, :T, :f, :ν, :ω])
    amplitude_info = test_field_parameters(ks, [:I, :E])

    c = 1u"c" |> u"m/s"

    @namespace!(field_params) do
        if !isempty([:λ, :T] ∩ ks)
            if :λ in ks
                T = λ/c
            elseif :T in ks
                λ = T/c
            end
            ν = 1/λ
            f = 1/T
            ω = 2π*u"rad"*f
        else # ∝ Frequency specified
            if :f in ks
                ν = f/c
                ω = 2π*u"rad"*f
            elseif :ν in ks
                f = ν*c
                ω = f/(2π*u"rad")
            elseif :ω in ks
                f = ω/(2π*u"rad")
                ν = f/c
            end
            T = 1/f
            λ = 1/ν
        end
    end()

    field_params
end

macro field(spec, var)
    spec.head == :-> ||
        error("Expected a block with parameters for definition of the field")
    block = spec.args[2]
    block.head == :block ||
        error("Expected a block with parameters for definition of the field")

    field_params = parse_block(block, Number)
    quote
        $(esc(var)) = make_field($field_params)
    end
end

export @field
