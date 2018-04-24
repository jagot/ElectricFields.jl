function parse_block(block, T)
    line_no = 0
    filename = ""
    error_message = v -> error("$(filename):$(line_no)\n>   $(v)")

    field_params = Dict{Symbol,T}()
    param_line_nos = Dict()

    for line in block.args
        typeof(line) == Expr || error_message("Expected expression, got $(line)")
        if line.head == :line
            line_no,filename = line.args
            continue
        end
        line.head == Symbol(":") ||
            error_message("Expected “parameter : value expression”, got $(line)")
        k = line.args[1]
        k in keys(field_params) &&
            error_message("Field parameter $(k) already specified at $(filename):$(param_line_nos[k])")
        v = line.args[2]
        field_params[k] = eval(v)
        param_line_nos[k] = line_no
    end

    field_params
end
