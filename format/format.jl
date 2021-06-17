using JuliaFormatter

function format_code()
    i = 0
    while i < 10
        if format(dirname(@__DIR__))
            return true
        end
        i += 1
    end
    return false
end

if format_code()
    @info("Formatting complete")
else
    @info("Formatting failed")
end
