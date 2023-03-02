# -*- coding: utf-8 -*-

from math import sqrt, radians, sin, cos, atan2, pi, asin

# Any Station instance, can have any (or all) of these attributes
station_member_names = ['name', 'lat', 'lon', 've', 'vn', 'se', 'sn', 'rho', 't']

class Station:
    '''A simple Station class.

        This module defines a Station class. A station is supposed to represent a
        point on the globe. It has coordinates (usually defined as longtitude and
        latitude), a name, (tectonic) velocities and respective standard deviations
        (in east and north components), a correlation coefficient between East and
        North velocity components and a time-span.
        This class is designed to assist the estimation of strain tensors; hence,
        only attributes that could help with this are considered.
        Этот модуль определяет класс Station. Станция должна представлять
         точка на глобусе. Он имеет координаты (обычно определяемые как долгота и
         широта), имя, (тектонические) скорости и соответствующие стандартные отклонения
         (по восточной и северной компонентам), коэффициент корреляции между востоком и
         Северные компоненты скорости и временной интервал.
         Этот класс предназначен для помощи в оценке тензоров деформации; следовательно,
         учитываются только те атрибуты, которые могли бы в этом помочь.

        Attributes:
            name (str) : the name of the station
            lon (float): longtitude of the station (radians). In case the station
                         coordinates are transformed to easting and northing
                         (aka to projection coordinates), this component will
                         hold the Easting.
            lat (float): latitude of the station (radians). In case the station
                         coordinates are transformed to easting and northing
                         (aka to projection coordinates), this component will
                         hold the Northing.

            ve (float) : velocity of the east component, in meters/year
            vn (float) : velocity of the north component, in meters/year
            se (float) : std. deviation of the east velocity component, in
                         meters/year
            sn (float) : std. deviation of the north velocity component, in
                         meters/year
            rho (float): correlation coefficient between East and North velocity
                         components
            t (float)  : time-span in decimal years
            
    '''

    def __init__(self, *args, **kargs):
        '''Station constructor.
        
            Station constructor; construction can be performed:
                #. from an input string of type:
                    "name lon lat Ve Vn Se Sn RHO T"
                    where lon and lat are in decimal degrees and velocity components
                    are in mm/year.
                #. given any of the (above mentioned) instance members/attributes.

            e.g. s = Station("akyr +24.91260690 +34.98083160 8.71244 -15.1236 0.00136367 0.000278371 0.5  2.5")
                 s = Station(name="akyr")
                 s = Station(name="akyr", lat=34.98083160, ve=-0.0151236)

            Args:
                *args (str): if provided, then it is supposed to be a station
                             string ("name lon lat Ve Vn Se Sn RHO T") and the
                             function will try to resolve it and assign member
                             values.
                **kargs:     any named member variable, aka one of:
                    * name
                    * lon
                    * lat
                    * ve
                    * vn
                    * se
                    * sn
                    * rho
                    * t
        '''
        self.set_none()

        if len(args) != 0: self.init_from_ascii_line(args[0])# если в конструкторе список не роавен 0 вызов init_from_ascii_line

        if len(kargs) != 0: # если именнованный словарик не равен 0
            for key, val in kargs.items():# распаковываем на ключ и значен
                if key in station_member_names:# если ключ в station_member_names = ['name', 'lat', 'lon', 've', 'vn', 'se', 'sn', 'rho', 't']
                   # print('vhod dannie ',key, val)
                    setattr(self, key, val)# устанавливает атрибуты ключ и значение

    def init_from_ascii_line(self, input_line):
        '''Assignment from string.

            This function will initialize all member values of a station
            instance, given a (string) line of type:
            "name lon lat Ve Vn Se Sn RHO T"
            where lon and lat are in decimal degrees and velocity components
            and sigmas are in mm/year.

            Args:
                input_line (str): a string (line) of type:
                                  "name lon lat Ve Vn Se Sn RHO T"

            Raises:
                RuntimeError: if the input line (string) cannot be resolved
        '''
        l = input_line.split()
        try:# Входныке данные переводит градусы в радианы скорости и ошибки делит на 1e3 т.е. +0.0
            self.name = l[0]
            self.lon  = radians(float(l[1]))
            self.lat  = radians(float(l[2]))
            self.ve   = float(l[3]) / 1e3
            self.vn   = float(l[4]) / 1e3
            self.se   = float(l[5]) / 1e3
            self.sn   = float(l[6]) / 1e3
            self.rho  = float(l[7]) / 1e3
            self.t    = float(l[8])
        except:
            print('[DEBUG] Invalid Station instance constrution.')
            print('[DEBUG] Input line \"{}\"'.format(input_line.strip()))
            raise RuntimeError

    def set_none(self):
        '''Set to None.

            Set all instance member values to None.
        '''
        self.name = None                                                        
        self.lon  = None
        self.lat  = None
        self.ve   = None
        self.vn   = None
        self.se   = None
        self.sn   = None
        self.rho  = None
        self.t    = None

    def distance_from(self, sta):
        '''Distance to another station.

            Compute the distance of the instance to another instance (of type
            Station). The algorithm is just an Eucledian norm, so the actual
            station components must already have been transformed to cartesian.

            Args:
                sta (Station): a station instance

            Returns:
                tuple (float, float, float): a 3-float tuple, where the elements
                are
                    #. dlon 
                    #. dlat
                    #. dr
                If the calling station has index i and the station passed in
                has index j, then the returned values are computed as
                    * δlon = lon_j - lon_i
                    * δlat = lat_j - lat_i
                    * δr   = sqrt{δlon**2 + δlat**2}
              
            Warning:
                For the function to return valid results, the station coordinate
                component must not be in ellipsoidal coordinates; the station
                coordinates should have already been transformed to cartesian
                before calling this function. The function will treat the "lon"
                attribute as "x" or "Easting" and the "lat" component as "y" or
                "Northing".

        '''
        print('station.py distance from')
        print('sta lon', sta.lon)
        print('self ', self.lon)
        print('sta lat ', sta.lat)
        print('self lat', self.lat)
        dlon = sta.lon - self.lon
        dlat = sta.lat - self.lat
        print('dlon', dlon)
        print('dlat',dlat)
        return dlon, dlat, sqrt(dlat*dlat + dlon*dlon)
    
    def squared_distance_from(self, sta):
        '''Squared distance to another station.

            Compute the squared distance of the instance to another instance (of type
            Station). The algorithm is just an Eucledian norm, so the actual
            station components must already have been transformed to cartesian.
            Note that the station coordinates (aka .lon and .lat) are first
            divided with 1e3, so that the squared distance is not too big a
            number.Расстояние в квадрате до другой станции.
             Вычислить квадрат расстояния экземпляра до другого экземпляра (типа
             Станция). Алгоритм представляет собой просто евклидову норму, поэтому фактический
             компоненты станции должны быть уже преобразованы в декартову систему.
             Обратите внимание, что координаты станции (иначе .lon и .lat) идут первыми.
             разделить на 1e3, чтобы квадрат расстояния не был слишком большим a
             количество.

            Args:
                sta (Station): a station instance

            Returns:
                    #. dr^2
                If the calling station has index i and the station passed in
                has index j, then the returned values are computed as
                Если вызывающая станция имеет индекс i и станция передана в
                 имеет индекс j, то возвращаемые значения вычисляются как
                    * δlon = (lon_j - lon_i)/1e3
                    * δlat = (lat_j - lat_i)/1e3
                    * δr   = δlon**2 + δlat**2
              
            Warning:
                For the function to return valid results, the station coordinate
                component must not be in ellipsoidal coordinates; the station
                coordinates should have already been transformed to cartesian
                before calling this function. The function will treat the "lon"
                attribute as "x" or "Easting" and the "lat" component as "y" or
                "Northing".

                This function is only useful for optimizing distance_from wrt
                time, as it does not call the sqrt function. Profile before you
                use and use with care!

        '''
        #print(sta.lon,'- ' ,self.name,' -  ', self.lon)
        dlon = (sta.lon - self.lon)/1e3
        dlat = (sta.lat - self.lat)/1e3
        return (dlat*dlat + dlon*dlon)

    def haversine_distance(self, sta, R=6372797.560856e0): 
#        self.d = 'ok'
#        
#        """Computes the distance, in meters, between two points on a sphere.
#
#            Determine the great-circle distance between two points on a sphere
#            given their longitudes and latitudes (self.lat, self.lon in radians)
#            see https://en.wikipedia.org/wiki/Haversine_formulaВычисляет расстояние в метрах между двумя точками на сфере.
#             Определить расстояние по дуге большого круга между двумя точками на сфере
#             учитывая их долготу и широту (self.lat, self.lon в радианах)
#           
#        """
#   
        def ArcInRadians(frm, to):
#            """Computes the arc, in radians, between two points on a sphere.
#               Вычисляет дугу в радианах между двумя точками на сфере
#            """
#            print('arcinRadians ',frm.lon)
            latarc = frm.lat - to.lat    #max -mean
            lonarc = frm.lon - to.lon    #max -mean
            lath   = sin(latarc*0.5e0)    #sin(max -mean)
            lath  *= lath                 #sin(max -mean)*sin(max -mean)  т.е.  sin(max -mean)^2
            lonh   = sin(lonarc*0.5e0)
            lonh  *= lonh
            tmp    = cos(frm.lat) * cos(to.lat)  #cos(max lat)* cos(mean Lat)
#            Vernet Rzemli*2*arcsin(sqrt(  sin(LATmax -LATmean)^2 +  (cos(max lat)* cos(mean Lat)) *  sin(LONmax -LONmean)^2  ))

            return 2e0 * asin(sqrt(lath + tmp*lonh))
#        return  sta.lat
#        return ArcInRadians(self, sta)
        return R*ArcInRadians(self, sta)
