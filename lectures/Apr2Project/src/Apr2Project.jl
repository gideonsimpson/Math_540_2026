module Apr2Project

# by default, this is not exported, so it's immediately visible
# if you want to use it, you have to do Apr2Project.greet()
greet() = print("Hello World!")

greet_person(name) = print("Hello $(name)!")
# exported, so you can use it without the module prefix
export greet_person

include("examplefunctions.jl")
export add, subtract
end # module Apr2Project
