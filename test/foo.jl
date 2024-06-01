
#=
writing a macro to add parameters to an existing paramstruct
=#

@defparams CoolParams :gamma 0.5

eval(genstruct(:MyParams2, :gamma, 0.5))

MyParams2().gamma


using DynamicHydra
p1 = DynamicHydraParams()
p2 = 

p = @with_kw (f = x -> x^3, y = 3, z = "foo")

p

@defparams NewParamStruct

@with_kw(a = 1, b = 3)


@edit @with_kw

NewParamStruct(1).Idot_max_rel


eval(xpr)
Params


macro extend(oldstruct::DataType, paramstructname::Symbol)
    quote
        mutable struct $paramstructname

        end
    end
end
 
begin
    macro extend_struct(struct_type, new_struct_name, fields...)
        quote
            @with_kw mutable struct $new_struct_name
                
                #$(Expr(:copyn, struct_type().parameters.args[2:end]...))
                #$(Expr(:tovector, Expr(:tuple, esc.(fields)...)))
            end
        end
    end
    
    @extend_struct DEBParams DEB2Params p2::Real = 2.
    
end



DEBParams(p1 = 2.)
