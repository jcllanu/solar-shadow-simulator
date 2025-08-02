import geometrical_utilities as geom
import math
import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os

EARTH_ROTATION_AXIS_ANGLE_DEG = 23.44
EARTH_ROTATION_AXIS_ANGLE_RAD = geom.degree2radian(23.44)
MONTHS = ['January', 'February', 'March', 'April', 'May', 'June',
            'July', 'August','September', 'October', 'November', 'December']
DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

def position(latitude, longitude, t):
    (x0, y0, z0) = geom.spherical2cartesian(1,latitude,longitude)
    (r0, alpha0) = geom.cartesian2polar(x0,y0)
    (x2, y2) = geom.polar2cartesian(r0, alpha0 + angle_rotation_earth(t))
    (x3, y3, z3) = geom.rot_fix_y(-EARTH_ROTATION_AXIS_ANGLE_RAD, x2, y2, z0)
    return (x3,y3,z3)

def angle_sun_earth(t): #t en días
    return (t-(datetime(2025, 6, 21) - datetime(2025, 1, 1)).days) *2*math.pi/365.25 # June 21st

def angle_rotation_earth(t): # t en días
    return 2*math.pi*t+ math.pi + (t-(datetime(2025, 6, 21) - datetime(2025, 1, 1)).days) *2*math.pi/365.25

def get_time_in_days(month, day, hour, minute):
    return sum(DAYS_IN_MONTH[:month-1]) + day + (hour + minute/60)/24

def init_function():
    return

def plot_shadow(month, day, hour, minute, latitude_rad, longitude_rad, ax, points, mode):
    ax.clear()
    t = get_time_in_days(month, day, hour, minute)
    
    normal_vector = np.array(position(latitude_rad, longitude_rad, t))

    angle_sun = angle_sun_earth(t)
    sun_ray = np.array((-math.cos(angle_sun), -math.sin(angle_sun), 0))
    rotation_axis_vector = np.array((math.sin(EARTH_ROTATION_AXIS_ANGLE_RAD), 0, math.cos(EARTH_ROTATION_AXIS_ANGLE_RAD)))

    north = np.dot(normal_vector,rotation_axis_vector)*normal_vector + rotation_axis_vector
    east = np.cross(rotation_axis_vector, normal_vector)

    dot_product = np.dot(normal_vector,sun_ray)
    if abs(dot_product)<1e-10:   
        # sunrise or sunset
        pass 
    elif dot_product > 0:
        # Night
        pass
    else:
        shadow_vector = normal_vector - sun_ray / dot_product
        shadow_east_component = float(np.dot(shadow_vector, east))
        shadow_north_component = float(np.dot(shadow_vector, north))
        points.append([shadow_east_component, shadow_north_component])
        if mode == 'DAY':
            for point in points:
                ax.plot([0, point[0]], [0, point[1]], linestyle='-', color='black')
        elif mode == 'HOUR':
            ax.plot([point[0] for point in points], [point[1] for point in points], linestyle='-', color='red')
        
        ax.plot([0, shadow_east_component], [0, shadow_north_component], linestyle='-', color='black')

        # Add a label
        (r, angle) = geom.cartesian2polar(shadow_east_component, shadow_north_component)
        label = 'x = ' + str(round(shadow_east_component, 2)) + ' y = ' + str(round(shadow_north_component, 2)) + '\n' +\
                        'r = ' + str(round(r,2)) + ' angle = '+ str(round(geom.radian2degree(angle),2)) + '°'
        ax.text(shadow_east_component + 1, shadow_north_component + 1, label, fontsize=12)
        ax.scatter(shadow_east_component, shadow_north_component, color='blue',s=20)
        ax.scatter(-shadow_east_component, -shadow_north_component, color='yellow', s=100)

    L = 4
    ax.set_ylim(-L, L)
    ax.set_xlim(-L, L)
    time_label = MONTHS[month-1] +' '+ str(day) + ' UTC '+str(hour).rjust(2,'0')+':'+str(minute).rjust(2,'0')
    
    ax.set_title(time_label)
    ax.set_xlabel('WEST - EAST')
    ax.set_ylabel('SOUTH - NORTH')
    position_label = 'Latitude: ' + str(round(geom.radian2degree(latitude_rad),2))+'° ' +\
        'Longitude: '+ str(round(geom.radian2degree(longitude_rad),2))+'°'
    plt.suptitle(position_label) 

def plot_shadow_day(time_0_1, month, day, latitude_rad, longitude_rad, ax, points):
    hours = time_0_1*24
    hour = math.floor(hours)
    minute = round((hours-hour)*60)
    plot_shadow(month, day, hour, minute, latitude_rad, longitude_rad, ax, points, 'DAY')

def plot_shadow_hour(day_year, hour, minute, latitude_rad, longitude_rad, ax, points):
    ax.cla()
    date = datetime(2025, 1, 1) + timedelta(days=day_year - 1)
    month, day = date.month, date.day
    plot_shadow(month, day, hour, minute, latitude_rad, longitude_rad, ax, points, 'HOUR')

def bisection_method(a, b, f, fa, fb, precision):

    if fa*fb >= 0:
        print('a = ' + str(a) + ' b = ' + str(b) + ' f(a) = ' + str(fa) + ' f(b) = ' + str(fb))
        return None
    
    c = (a+b)/2
    fc = f(c)
    if abs(fc) < precision:
        return c
    elif fc > 0:
        if fa > fb:
            return bisection_method(c, b, f, fc, fb, precision)
        else:
            return bisection_method(a, c, f, fa, fc, precision)

    else:
        if fa > fb:
            return bisection_method(a, c, f, fa, fc, precision)
        else:
            return bisection_method(c, b, f, fc, fb, precision)
            
def get_sunrise_and_sunset_time(latitude_rad, month, day):
    year_day = get_time_in_days(month, day, 0, 0)

    f = lambda t:  np.dot(np.array(position(latitude_rad, 0, year_day + t)),
                            np.array((-math.cos(angle_sun_earth(year_day + t)), -math.sin(angle_sun_earth(year_day + t)), 0)))
    
    midnight = 0
    noon = 0.5 
    midnight_next_day = 1

    sunrise = bisection_method(midnight, noon, f, f(midnight), f(noon), 1/(24*60))
    sunset = bisection_method(noon, midnight_next_day, f, f(noon), f(midnight_next_day), 1/(24*60))

    hours = sunrise*24
    hour = math.floor(hours)
    minute = round((hours-hour)*60)
    print('\nLatitude: ' + str(round(geom.radian2degree(latitude_rad),2)) + '° ')
    print(MONTHS[month-1] +' '+ str(day))
    print('Sunrise: '+str(hour).rjust(2,'0')+':'+str(minute).rjust(2,'0')+ ' UTC')
    hours = sunset*24
    hour = math.floor(hours)
    minute = round((hours-hour)*60)
    print('Sunset: '+str(hour).rjust(2,'0')+':'+str(minute).rjust(2,'0')+ ' UTC')

def plot_sunrise_and_sunset_times(latitude_rad):
    sunrise=[]
    sunset=[]
    for day in range(365):

        f = lambda t:  np.dot(np.array(position(latitude_rad, 0, day + t)),
                               np.array((-math.cos(angle_sun_earth(day + t)), -math.sin(angle_sun_earth(day + t)), 0)))
        
        midnight = 0
        noon = 0.5 
        midnight_next_day = 1
        sunrise.append(bisection_method(midnight, noon, f, f(midnight), f(noon), 1/(24*60)))
        sunset.append(bisection_method(noon, midnight_next_day, f, f(noon), f(midnight_next_day), 1/(24*60)))
    
    days = 365
    start_date = datetime(2025, 1, 1)
    dates = [start_date + timedelta(days=i) for i in range(days)]

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot lines
    ax.plot(dates, sunset, label='Sunset', color='purple', linewidth=2)
    ax.plot(dates, sunrise, label='Sunrise', color='orange', linewidth=2)
    # Add vertical lines for solstices and equinoxes
    event_dates = {
        "Spring Equinox": datetime(2025, 3, 20),
        "Summer Solstice": datetime(2025, 6, 21),
        "Autumn Equinox": datetime(2025, 9, 22),
        "Winter Solstice": datetime(2025, 12, 21)
    }

    for label, date in event_dates.items():
        ax.axvline(date, color='gray', linestyle='--', linewidth=1)
        ax.text(date + timedelta(days=-5), 0.4, label, rotation=90, va='bottom', ha='center', fontsize=9, color='gray')
    

    # Format x-axis
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

    # Format y-axis to show time labels
    yticks = np.linspace(0, 1, 25)  # 0 to 24h in 3-hour steps
    ytick_labels = [f"{int(hours):02d}:{int(0):02d}" for hours in range(25)]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels)

    # Titles and labels
    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Time of Day (UCT)', fontsize=12)
    plt.suptitle('Sunrise and Sunset Times Over the Year\n\n'+'Latitude: ' + str(round(geom.radian2degree(latitude_rad),2)) + '°') 

    # Aesthetics
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend()
    fig.tight_layout()

    plt.show()
            
def generate_shadow_gif_given_hour(latitude_deg, longitude_deg, hour, minute):
    figure, ax = plt.subplots()
    latitude_rad = geom.degree2radian(latitude_deg)
    longitude_rad = geom.degree2radian(longitude_deg)
    

    print('Latitude: ' + str(round(latitude_deg)) + '° Longitude: ' + str(round(longitude_deg)) + '° ' + f"{int(hour):02d}:{int(minute):02d}")
    points=[]
    ani = animation.FuncAnimation(
            figure,
            func = plot_shadow_hour,
            frames=[i for i in range(366)],
            fargs=(hour, minute, latitude_rad, longitude_rad, ax, points),
            init_func=init_function
        )

    ani.save('lat'+str(round(latitude_deg))+'_long'+ str(round(longitude_deg))+ f"_{int(hour):02d}_{int(minute):02d}" +'.gif', writer='pillow', fps=100) 

def generate_shadow_gif_given_day(latitude_deg, longitude_deg, month, day):
    figure, ax = plt.subplots()
    latitude_rad = geom.degree2radian(latitude_deg)
    longitude_rad = geom.degree2radian(longitude_deg)
    

    print('Latitude: ' + str(round(latitude_deg)) + '° Longitude: ' + str(round(longitude_deg)) + '° ' + MONTHS[month-1] +' ' + str(day))
    points=[]
    ani = animation.FuncAnimation(
            figure,
            func = plot_shadow_day,
            frames=np.linspace(0, 1, 24*6*2 + 1),
            fargs=(month, day, latitude_rad, longitude_rad, ax, points),
            init_func=init_function
        )

    ani.save('lat'+str(round(latitude_deg))+'_long'+ str(round(longitude_deg))+MONTHS[month-1] + str(day)+'.gif', writer='pillow', fps=100) 

def generate_shadow_day_gifs_all_latitudes():
    figure, ax = plt.subplots()
    dates = [(2,3), (3,20), (5,4), (6,21), (8,6), (9,22), (11,6), (12,21)]
    for i in range(1,18):
        count = 0
        for date in dates:
            count += 1
            latitude_deg = 90 - 10 * i
            longitude_deg = 0
            latitude_rad = geom.degree2radian(latitude_deg)
            longitude_rad = geom.degree2radian(longitude_deg)
            month = date[0]
            day = date[1]
            print('Latitude: ' + str(round(latitude_deg)) + '° Longitude: ' + str(round(longitude_deg)) + '° ' + MONTHS[month-1] +' ' + str(day))

            points=[]
            ani = animation.FuncAnimation(
                    figure,
                    func = plot_shadow_day,
                    frames=np.linspace(0, 1, 24*6*2 + 1),
                    fargs=(month, day, latitude_rad, longitude_rad, ax, points),
                    init_func=init_function
                )
            directory = 'shadow_evolution/lat'+str(round(latitude_deg))+'/fixed_day/'
            os.makedirs(directory, exist_ok=True)  # Create folder if it doesn't exist
            ani.save(directory + str(count)+'_lat'+str(round(latitude_deg))+'_long'+ str(round(longitude_deg))+MONTHS[month-1] + str(day)+'.gif', writer='pillow', fps=100) 

def generate_shadow_hour_gifs_all_latitudes():
    figure, ax = plt.subplots()
    minute = 0
    for i in range(1,18):
        count = 0
        for hour in [6, 8, 10, 12, 14, 16, 18]:
            count += 1
            latitude_deg = 90 - 10 * i
            longitude_deg = 0
            latitude_rad = geom.degree2radian(latitude_deg)
            longitude_rad = geom.degree2radian(longitude_deg)
            print('Latitude: ' + str(round(latitude_deg)) + '° Longitude: ' + str(round(longitude_deg)) + '° ' + f"{int(hour):02d}:{int(minute):02d}")
            points=[]
            ani = animation.FuncAnimation(
                    figure,
                    func = plot_shadow_hour,
                    frames=[i for i in range(366)],
                    fargs=(hour, minute, latitude_rad, longitude_rad, ax, points),
                    init_func=init_function
                )
            
            directory = 'shadow_evolution/lat'+str(round(latitude_deg))+'/fixed_hour/'
            os.makedirs(directory, exist_ok=True)  # Create folder if it doesn't exist
            ani.save(directory + str(count)+'_lat'+str(round(latitude_deg))+'_long'+ str(round(longitude_deg))+ f"_{int(hour):02d}_{int(minute):02d}" +'.gif', writer='pillow', fps=100) 

            
if __name__=="__main__":
    mode = int(input('1.- Generate sun trajectory and shadow evolution during a daylength in informative latitudes for solstices, equinoxes, and specific days in between \n' +\
                 '2.- Generate sun trajectory and shadow evolution during a given day in a specific location\n'+\
                 '3.- Generate sun and shadow evolution throughout a year in a set of informative hours and latitudes \n'+\
                 '4.- Generate sun and shadow evolution throughout a year in a given time of the day and a specific location\n'+\
                 '5.- Plot sunrise and sunset times during the year for a given latitude\n'+\
                 '6.- Plot sunrise and sunset times of a given day for a given latitude\n'+\
                 'Select funcionality: '))
    if mode == 1:
        generate_shadow_day_gifs_all_latitudes()
    elif mode == 2:
        latitude_deg = float(input('Insert latitude (in degrees) as a decimal number (e.g. 40.4170): '))
        longitude_deg = float(input('Insert longitude (in degrees) as a decimal number (e.g. -3.7034): '))
        month = int(input('Insert the month as a number (e.g. January - 1, February - 2 ...): '))
        day = int(input('Insert the day as a number: '))
        generate_shadow_gif_given_day(latitude_deg, longitude_deg, month, day)
    elif mode == 3:
        generate_shadow_hour_gifs_all_latitudes()
    elif mode == 4:
        latitude_deg = float(input('Insert latitude (in degrees) as a decimal number (e.g. 40.4170): '))
        longitude_deg = float(input('Insert longitude (in degrees) as a decimal number (e.g. -3.7034): '))
        hour = int(input('Insert the hour of the day as an integer from 0 to 23: '))
        minute = int(input('Insert the minute in that hour as an integer from 0 to 59: '))
        generate_shadow_gif_given_hour(latitude_deg, longitude_deg, hour, minute)
    elif mode == 5:
        latitude_deg = float(input('Insert latitude (in degrees) as a decimal number (e.g. 40.4170): '))
        plot_sunrise_and_sunset_times(geom.degree2radian(latitude_deg))
    elif mode == 6:
        latitude_deg = float(input('Insert latitude (in degrees) as a decimal number (e.g. 40.4170): '))
        month = int(input('Insert the month as a number (e.g. January - 1, February - 2 ...): '))
        day = int(input('Insert the day as a number: '))
        get_sunrise_and_sunset_time(geom.degree2radian(latitude_deg), month, day)
        