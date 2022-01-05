 function thisFun = func1(x, y, e0)

    f1 = y(2);
    f2 = (x^2 -e0 ) * y(1);

    thisFun = [f1; f2];
    end
