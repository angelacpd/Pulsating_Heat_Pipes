function [x_out, y_out] = get_local_max(Y)
    
    %window = 10;
    len_v = length(Y);
    x = 1:len_v;
    
    TF = islocalmax(Y)
    x_out = (x(TF))'
    y_out = Y(TF)
    plot(x,Y,x(TF),Y(TF),'r*')
    
end