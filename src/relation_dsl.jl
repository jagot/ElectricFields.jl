# * Parse assignment DSL blocks

function parse_block(block, T)
    line_no = 0
    filename = ""
    error_message = v -> error("$(filename):$(line_no)\n>   $(v)")

    field_params = Dict{Symbol,T}()
    param_line_nos = Dict()

    for line in block.args
        if typeof(line) == LineNumberNode
            line_no,filename = line.line,line.file
            continue
        end
        typeof(line) == Expr || error_message("Expected expression, got $(line)")
        line.head == Symbol("=") ||
            error_message("Expected “parameter = value expression”, got $(line)")
        k = line.args[1]
        k in keys(field_params) &&
            error_message("Field parameter $(k) already specified at $(filename):$(param_line_nos[k])")
        v = line.args[2]
        field_params[k] = eval(v)
        param_line_nos[k] = line_no
    end

    field_params
end

# * DSL for calculation of quantities

# The somewhat complicated setup with walking the expression tree in
# =make_field= (twice implemented :) allows for very clean expression
# of the formulaic dependencies between different quantities, almost
# as if it were pure math. All ingoing quantities are either unitful,
# or made unitful using the set base units. Furthermore, outgoing
# quantities are transformed to the base units, even though the
# expression may result in a different (but equivalent) unit
# expression. This way, even if a period time is provided in
# femtoseconds, the wavenumber will always be returned in Kaysers, for
# instance.

# ** Testing presence of quantity in the namespace
# These methods are used for the =if x= constructs in the quantity
# conversion DSL. We return =true= if the symbol =x= is present in
# the =params= dict. Since we implement it as an expression walking
# algorithm, we can have constructs as =if x || y=, which will
# expand to =if :x in keys(params) || :y in keys(params)=. At first,
# however, we check if the symbols =isdefined= in the =Main= scope,
# which is generally the case for operators and such, to avoid
# generating tentative key lookups for every symbol that is
# encountered.

test_symbol_walk(node, params) = node
test_symbol_walk(node::Symbol, params) =
    isdefined(Main, node) ? node : Expr(:call, :in, Expr(:quote, node),
                                        Expr(:call, :keys, params))
test_symbol_walk(node::Expr, params) =
    Expr(node.head, test_symbol_walk.(node.args, Ref(params))...)# 

# ** Expression walker
# This walks the quantity conversions block, replacing symbols with
# references to dictionary items and converting quantities to
# applicable base units, before assignment. As above, before
# generating dictionary lookups for all symbols encountered, we
# first check if the symbol =isdefined= in the =Main= scope.

get_reference(r::Symbol,params) = Expr(:ref, params, Expr(:quote, r))
get_reference(r,params) = r

function get_symbol(s::Symbol, params, escape)
    if isdefined(Main, s)
        escape ? esc(s) : s
    else
        test_expression = test_symbol_walk(s, params)
        ref_expression = get_reference(s, params)
        Expr(:if, test_expression, ref_expression, escape ? esc(s) : s)
    end
end
get_symbol(s, params, escape) = s

walk(node, params, escape) = node
walk(node::Symbol, params, escape) = get_symbol(node, params, escape)

function walk(node::Expr, params, escape)
    if node.head ∈ [:line, :quote, :macrocall]
        node
    elseif node.head ∈ [:if, :elseif]
        # Dispatch "if x" to test_symbol_walk, which checks if :x is
        # present as key in params.
        args = Any[test_symbol_walk(node.args[1], params)]
        append!(args, [walk(a, params, escape) for a in node.args[2:end]])
        Expr(node.head, args...)
    else
        if node.head == :(=)
            target,expr = node.args
            r = get_reference(target, params)
            result = walk(expr, params, escape)
            Expr(node.head, r, result)
        else
            args = walk.(node.args, Ref(params), Ref(escape))
            Expr(node.head, args...)
        end
    end
end

# ** Namespace macro

"""
    @namespace!(exprs, params)

This macro uses the dictionary `params` as a "namespace", i.e. all
symbols are assumed to be keys in this dictionary. We need to escape
the generated expression tree to modify in the scope of the caller,
not the global scope.
"""
macro namespace!(exprs, params)
    local tree = walk(exprs, params, false)
    quote
        $(esc(tree))()
    end
end

# ** Test of "competing quantities"
"""
    test_field_parameters(field_params, set)

This function ensures that one and only one of \"competing\"
quantities is specified."""
function test_field_parameters(field_params, set)
    info = set ∩ keys(field_params)
    set_string = join(set, ", ", " and ")

    length(info) == 0 &&
        (length(set) > 1 && error("Need to provide one of $(set_string)") ||
         error("Required parameter $(set_string) missing"))
    length(info) > 1 &&
        error("Can only specify one of $(set_string)")

    info
end
