classdef wave;
    % Wavefront of a phase with given parameters
    properties
        lambda;
        height_error;
        amp;
        wavefront_error;
        front;
    end
    
    methods
        function obj = wave(lambda,height_error, amp,g)
            if (nargin > 0)
                obj.lambda = lambda;
                obj.height_error = height_error;
                obj.amp = amp;
                obj.wavefront_error = 4*pi*height_error/lambda;
                obj.front = g.circa.*amp*exp(1i*obj.wavefront_error);
            end
        end
        
        function obj = update(obj,g)
            obj.wavefront_error = 4*pi*obj.height_error/obj.lambda;
            obj.front = g.circa.*obj.amp.*exp(1i*obj.wavefront_error);
        end
        
        function obj = aberration(obj,n,g)
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
            obj.height_error = eval(n);
            obj = obj.update(g);
        
        end
        function obj = plot_height(obj,D,g)
            if (D==2)
                xh = -1:g.resolution:1;
                image(xh,xh,g.circa.*obj.height_error,'CDataMapping','scaled')
                set(gca,'visible','off')
                colormap(gray)
                shading interp
            else
                surfl(g.x,g.y,g.circa.*obj.height_error)
                set(gca,'visible','off')
                colormap(gray)
                shading interp 
            end
        end
        
        function obj = tilt(obj,x_mag,y_mag,g)
            obj.height_error = (x_mag.*g.x + y_mag.*g.y);
            obj = obj.update(g);
        end
    end
end
