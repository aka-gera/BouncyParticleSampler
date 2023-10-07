function ff=Reflect(gradient, v)

    ff= v - 2 * ( (gradient)' * v / dot(gradient,gradient)) * gradient;

end
