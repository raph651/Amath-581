 function thisFun = func3(x, y, e0, g)

    f1 = y(2);
    f2 = (g*abs(y(1))^2 +x^2-e0 ) * y(1);

    thisFun = [f1; f2];
    end
