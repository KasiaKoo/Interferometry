classdef height_error
    properties
        distribution;
        shape;
        rms;
    end
    methods
        function obj = height_error(g,func)
            if (nargin > 0)  
                if func(1) == 'Z'
                    obj = obj.Zernike(func(2),g);
                elseif func(1) == 'T'
                    obj = obj.tilt(func(2),func(3),g);
                elseif func(1) == 'F'
                    obj.distribution = 1;
                elseif func(1) == 'R'
                    obj = obj.random(g);
                    obj.rms = sqrt(nanmean(obj.distribution.*g.circa(:).^2) ...
                    -nanmean(obj.distribution.*g.circa(:)));
                end
                
            end 
        end
                
        function obj = Zernike(obj,n,g)
            Z_00 = 1;
            Z_11 = 2.*g.r.*cos(g.theta);
            Z_1n1 = 2.*g.r.*sin(g.theta);
            Z_20 = sqrt(3).*(2.*(g.r.^2) - 1);
            Z_22 = sqrt(6).*(g.r.^2).*cos(2.*g.theta);
            Z_2n2 = sqrt(6).*(g.r.^2).*sin(2.*g.theta);
            Z_31 = 2*sqrt(2).*(3.*(g.r.^3) - 2.*g.r).*cos(g.theta);
            Z_3n1 = 2*sqrt(2).*(3.*(g.r.^3) - 2.*g.r).*sin(g.theta);
            Z_33 = 2*sqrt(2).*(g.r.^3).*cos(3.*g.theta);
            Z_3n3 = 2*sqrt(2).*(g.r.^3).*sin(3.*g.theta);
            Z_40 = sqrt(5).*(6.*(g.r.^4) - 6.*(g.r.^2)+1);
            obj.distribution = eval(n);
        end
        
        function obj = tilt(obj,x_mag,y_mag,g)
            obj.distribution = (x_mag.*g.x + y_mag.*g.y);
        end
        
        function obj = random(obj,g)
            ri = rand( [50 50])-0.5;
            kern = imresize(exp(-g.r.^2/(2*0.25^2)),[30 30]);
            ri = conv2(ri,kern);
            ri = imresize(ri,20,'bicubic');
            obj.distribution = ri(300:700,300:700);
        end
        
        function obj = plot_height(obj,g)
            circular_dist = obj.distribution.*g.circa;
            figure();
            surfl(g.x,g.y,circular_dist);
            set(gca,'visible','off');
            colormap(gray);
            shading interp ;
            figure();
            imshow(circular_dist,[]);
        end
        
        function obj = zernike_coefficients(obj,g)
            Z_00 = 1.*g.circa;
            Z_11 = 2.*g.r.*cos(g.theta).*g.circa;
            Z_1n1 = 2.*g.r.*sin(g.theta).*g.circa;
            Z_20 = sqrt(3).*(2.*(g.r.^2) - 1).*g.circa;
            Z_22 = sqrt(6).*(g.r.^2).*cos(2.*g.theta).*g.circa;
            Z_2n2 = sqrt(6).*(g.r.^2).*sin(2.*g.theta).*g.circa;
            Z_31 = 2*sqrt(2).*(3.*(g.r.^3) - 2.*g.r).*cos(g.theta).*g.circa;
            Z_3n1 = 2*sqrt(2).*(3.*(g.r.^3) - 2.*g.r).*sin(g.theta).*g.circa;
            Z_33 = 2*sqrt(2).*(g.r.^3).*cos(3.*g.theta).*g.circa;
            Z_3n3 = 2*sqrt(2).*(g.r.^3).*sin(3.*g.theta).*g.circa;
            Z_40 = sqrt(5).*(6.*(g.r.^4) - 6.*(g.r.^2)+1).*g.circa;
            circular_dist = obj.distribution.*g.circa;
            c_nm = linspace(0,0,11);
            abb_options = {Z_00, Z_11, Z_1n1, Z_20, Z_22, Z_2n2, Z_31, Z_3n1, Z_33, Z_3n3, Z_40};
            for i=1:11
                abe = abb_options(i);
                c_nm(i) = sum(circular_dist(:).*abe(:))/sum(g.circa(:));
            end
            obj.c_nm = c_nm;
        end
        
  
    end
end

        