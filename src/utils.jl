export print_mem_address
export lrange

function print_mem_address(a)
    repr(UInt64(pointer_from_objref(a)))
end

function lrange(from::Float64, to::Float64, length::Int64 = 6)
    exp.(collect(range(log(from), log(to); length = length)))
end


