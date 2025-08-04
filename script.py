import geometrical_utilities as geom
import math
import numpy as np
from matplotlib import animation
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature

EARTH_ROTATION_AXIS_ANGLE_DEG = 23.44
EARTH_ROTATION_AXIS_ANGLE_RAD = geom.degree2radian(23.44)
MONTHS = ['January', 'February', 'March', 'April', 'May', 'June',
            'July', 'August','September', 'October', 'November', 'December']
DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

def position(latitude, longitude, t):
    """
    Calculate the 3D Cartesian coordinates of a point on Earth's surface at time t,
    taking into account Earth's rotation around its tilted axis.

    Parameters:
    - latitude (radians): geographic latitude of the point
    - longitude (radians): geographic longitude of the point
    - t (float): time in days, used to compute Earth's rotation angle

    Returns:
    - (x, y, z): tuple of Cartesian coordinates representing the position in 3D space
    """
    # Transform spherical coordinates of latitude and longitud to cartesian coordinates
    (x0, y0, z0) = geom.spherical2cartesian(1,latitude,longitude)  
    # Consider the plane z=z0 and get the polar coordinates of the point (x0,y0)
    (r0, alpha0) = geom.cartesian2polar(x0,y0)
    # Rotate the point (x0, y0) over the Z axis an angle which depends on the time t
    (x2, y2) = geom.polar2cartesian(r0, alpha0 + angle_rotation_earth(t))
    #Rotate the resulting point over the Y axis to account for the tilt of Earth's rotation axis
    (x3, y3, z3) = geom.rot_fix_y(-EARTH_ROTATION_AXIS_ANGLE_RAD, x2, y2, z0)
    return (x3,y3,z3)

def angle_sun_earth(t):
    """
    Calculate the angular position of the Earth around the Sun in radians at day t.

    Parameters:
    - t (float): time in days since the start of the year (e.g., January 1 = 0)

    Returns:
    - float: Earth's orbital angle in radians, where 0 corresponds to the reference date (June 21st, summer solstice)
    """
    return (t-(datetime(2025, 6, 21) - datetime(2025, 1, 1)).days) *2*math.pi/365

def angle_rotation_earth(t): 
    """
    Calculate the Earth's rotation angle around its own axis in radians at time t,
    including the effect of Earth's orbital position around the Sun.

    Parameters:
    - t (float): time in days since the start of the year (e.g., January 1 = 0)

    Returns:
    - float: Earth's rotation angle in radians, combining daily rotation and seasonal orbital adjustment.
    """
    return 2*math.pi*t+ math.pi + (t-(datetime(2025, 6, 21) - datetime(2025, 1, 1)).days) *2*math.pi/365

def location_in_which_sunrays_orthogonal(t):
    """
    Calculate the geographic location (latitude, longitude) on Earth where the Sun's rays are orthogonal
    at a given time t

    Parameters:
        t: Time parameter in days
    
    Returns:
        A tuple containing:
            - Latitude (in degrees) where the sun rays are orthogonal
            - Longitude (in degrees) where the sun rays are orthogonal
    """
    # Compute the angle between the Earth and the Sun. This angle is defined such that:
    #   - At the summer solstice (June 21), the angle is 0 radians.
    #   - At the autumn equinox, the angle is π/2 radians (90 degrees).
    angle_sun = angle_sun_earth(t)

    # Given the Earth-Sun angle, construct the unit vector representing sunlight direction,
    # but in the direction Earth to Sun. The sun ray lies in the XY-plane and is normal to Earth's surface
    (x0, y0, z0) = np.array((math.cos(angle_sun), math.sin(angle_sun), 0))

    #Rotate the resulting point over the Y axis to account for the tilt of Earth's rotation axis
    (x1, y1, z1) = geom.rot_fix_y(EARTH_ROTATION_AXIS_ANGLE_RAD, x0, y0, z0)
    
    # Consider the plane z=z1 and get the polar coordinates of the point (x1,y1)
    (r0, alpha0) = geom.cartesian2polar(x1,y1)

    # Rotate the point (x1, y1) over the Z axis an angle which depends on the time t
    (x2, y2) = geom.polar2cartesian(r0, alpha0 - angle_rotation_earth(t))

    # Transform cartesian coordinates into spherical coordinates (latitude and longitud)
    (r, latitude_rad, longitude_rad) = geom.cartesian2spherical(x2, y2, z1)
    return (geom.radian2degree(latitude_rad), geom.radian2degree(longitude_rad))

def get_time_in_days(month, day, hour, minute):
    """
    Convert a given date and time into a fractional day count within the year.

    Parameters:
    - month (int): month number (1-12)
    - day (int): day of the month
    - hour (int): hour of the day (0-23)
    - minute (int): minute of the hour (0-59)

    Returns:
    - float: time expressed as the number of days (including fractional days) since the start of the year
    """
    return sum(DAYS_IN_MONTH[:month-1]) + day + (hour + minute/60)/24

def bisection_method(a, b, f, fa, fb, precision):
    """
    Find a root of the function f in the interval [a, b] using the bisection method.

    Parameters:
    - a, b: Interval endpoints where f(a) and f(b) have opposite signs (root bracket).
    - f: Function for which to find the root.
    - fa, fb: Precomputed values f(a) and f(b).
    - precision: Desired precision for the root (stopping criterion based on |f(c)|).

    Returns:
    - c: Approximate root of f such that |f(c)| < precision.
    - None: If f(a) and f(b) do not bracket a root (i.e., have the same sign).

    Notes:
    - Recursively bisects the interval, choosing subintervals that contain the root
      based on the sign of function values at endpoints.
    - Stops when the function value at midpoint c is within the desired precision.
    """
    if fa*fb >= 0:
        print('a = ' + str(a) + ' b = ' + str(b) + ' f(a) = ' + str(fa) + ' f(b) = ' + str(fb))
        return None
    
    c = (a+b)/2 # midpoint
    fc = f(c) # precompute value: only one evaluation
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
   
def init_function():
    return

def plot_shadow(month, day, hour, minute, latitude_rad, longitude_rad, ax, points, mode):
    """
    Plot the shadow direction on a 2D plane for a given date, time, and geographic location.

    Parameters:
    - month (int): Month of the year (1-12)
    - day (int): Day of the month
    - hour (int): Hour of the day (0-23)
    - minute (int): Minute of the hour (0-59)
    - latitude_rad (float): Latitude in radians
    - longitude_rad (float): Longitude in radians
    - ax (matplotlib.axes.Axes): Matplotlib axis to plot on
    - points (list): List to accumulate shadow points for plotting paths
    - mode (str): Mode of plotting, 'DAY' plots individual shadow lines,
                  'HOUR' plots continuous path of shadows over time

    Behavior:
    - Clears the plot axis and calculates Earth's rotation and sun position
    - Computes the shadow vector based on sun rays and surface normal
    - Adds shadow direction vectors and labels to the plot
    - Handles edge cases like sunrise/sunset and nighttime (no shadow)
    - Updates plot limits, titles, and labels accordingly
    """
    ax.clear() # Clears the plot axis
    t = get_time_in_days(month, day, hour, minute) 

    # Calculate the position of the point defined by the latitude and longitude at time t.
    # This position lies on the unit sphere and can be interpreted as the surface normal vector
    # to the Earth at the given location and time.
    normal_vector = np.array(position(latitude_rad, longitude_rad, t)) 

    # Compute the angle between the Earth and the Sun. This angle is defined such that:
    #   - At the summer solstice (June 21), the angle is 0 radians.
    #   - At the autumn equinox, the angle is π/2 radians (90 degrees).
    angle_sun = angle_sun_earth(t)

    # Given the Earth-Sun angle, construct the vector representing sunlight direction 
    # from the Sun toward the Earth. The sun ray lies in the XY-plane.
    sun_ray = np.array((-math.cos(angle_sun), -math.sin(angle_sun), 0))

    # Define the Earth's rotation axis vector in 3D space.
    # It is tilted by approximately 23.44° from the Z-axis.
    rotation_axis_vector = np.array((
        math.sin(EARTH_ROTATION_AXIS_ANGLE_RAD), 
        0, 
        math.cos(EARTH_ROTATION_AXIS_ANGLE_RAD)
    ))

    # Define the local geographic directions (north and east):

    # --- NORTH VECTOR ---
    # The north vector lies at the intersection of:
    #   1. The tangent plane to the Earth's surface at the location:
    #        (This plane is orthogonal to normal_vector) → normal_vector · x = 0
    #   2. The plane spanned by the rotation axis and normal_vector:
    #        x = λ * normal_vector + μ * rotation_axis_vector
    #
    # To ensure the north vector points in the upward hemisphere (positive rotation_axis component),
    # we fix μ = 1 and compute the vector. It is then normalized to unit length.
    north_direction = -np.dot(normal_vector, rotation_axis_vector) * normal_vector + rotation_axis_vector
    north = north_direction / np.linalg.norm(north_direction)

    # --- EAST VECTOR ---
    # The east vector lies in:
    #   1. The tangent plane (normal_vector · x = 0)
    #   2. The plane perpendicular to the rotation axis (rotation_axis_vector · x = 0)
    #
    # This vector is obtained by taking the cross product of the rotation axis and normal vector.
    # By the right-hand rule, the resulting vector points east.
    east = np.cross(rotation_axis_vector, normal_vector)

    # Project the shadow if sunlight reaches the location:
    dot_product = np.dot(normal_vector, sun_ray)

    if abs(dot_product) < 1e-10:
        # The Sun is on the horizon — either sunrise or sunset
        pass 
    elif dot_product > 0:
        # The location is in Earth's shadow (nighttime)
        pass
    else:
        # The location is in daylight — compute the shadow direction
        # The shadow vector is computed using projection theory: 
        # Assume a unit stick orthogonal to the tangent plane, i.e., with its tip in normal_vector.
        # sun_ray (d) is the direction of incoming sunlight.
        # The shadow is where the line starting at normal_vector in direction sun_ray
        # intersects the tangent plane perpendicular to normal_vector.
        #
        # Line: p(t) = normal_vector + t * sun_ray
        # Tangent plane condition: normal_vector · p(t) = 0
        # Solve for t: t = -dot(normal_vector, normal_vector) / dot(normal_vector, sun_ray) = -1 / dot_product
        #
        # Shadow vector = p(t) = normal_vector - (1 / dot_product) * sun_ray
        shadow_vector = normal_vector - (1 / dot_product) * sun_ray

        # Project the shadow vector onto the east and north axes (both are unit vectors)
        shadow_east_component = float(np.dot(shadow_vector, east))
        shadow_north_component = float(np.dot(shadow_vector, north))

        # Add the new shadow point to the points list
        points.append([shadow_east_component, shadow_north_component])
        if mode == 'DAY': # Plot the evolution of the stick's shadow for the previous times
            for point in points:
                ax.plot([0, point[0]], [0, point[1]], linestyle='-', color='black')
        elif mode == 'HOUR': # Plot the evolution of the stick's tip shadow position for previous times
            ax.plot([point[0] for point in points], [point[1] for point in points], linestyle='-', color='red')
        
        # Plot the shadow for the current time
        ax.plot([0, shadow_east_component], [0, shadow_north_component], linestyle='-', color='black')

        # Transform the  cartesian coordinates in the east-north plane to polar coordinates
        (r, angle) = geom.cartesian2polar(shadow_east_component, shadow_north_component)
        # Add a label
        label = 'x = ' + str(round(shadow_east_component, 2)) + ' y = ' + str(round(shadow_north_component, 2)) + '\n' +\
                        'r = ' + str(round(r,2)) + ' angle = '+ str(round(geom.radian2degree(angle),2)) + '°'
        ax.text(shadow_east_component + 1, shadow_north_component + 1, label, fontsize=12)
        # Plot origin and stick's tip shadow position for current t
        ax.scatter(shadow_east_component, shadow_north_component, color='blue',s=20)
        ax.scatter(0, 0, color='blue',s=20)

        # Plot sun position at current time
        ax.scatter(-shadow_east_component, -shadow_north_component, color='yellow', s=100)
    
    # Set axis labels, limits and titles
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
    """
    Convert fractional day time to hour and minute, then plot the shadow for daytime.

    Parameters:
    - time_0_1: Time of day as a fraction between 0 (midnight) and 1 (next midnight).
    - month (int): month number (1-12)
    - day (int): day of the month
    - latitude_rad: Latitude in radians.
    - longitude_rad: Longitude in radians.
    - ax: Matplotlib axis to plot on.
    - points (list): List to accumulate shadow points for plotting paths

    The function converts the fractional time into discrete hour and minute values,
    then calls the `plot_shadow` function with the 'DAY' mode to plot shadows 
    corresponding to that time and location.
    """
    hours = time_0_1 * 24
    hour = math.floor(hours)
    minute = round((hours-hour)*60)
    plot_shadow(month, day, hour, minute, latitude_rad, longitude_rad, ax, points, 'DAY')

def plot_shadow_hour(day_year, hour, minute, latitude_rad, longitude_rad, ax, points):
    """
    Convert day of year to month and day, then plot the shadow for a specific hour and minute.

    Parameters:
    - hour (int): Hour of the day (0-23)
    - minute (int): Minute of the hour (0-59)
    - latitude_rad (float): Latitude in radians
    - longitude_rad (float): Longitude in radians
    - ax (matplotlib.axes.Axes): Matplotlib axis to plot on
    - points (list): List to accumulate shadow points for plotting paths


    The function converts the day of the year to a calendar date (month and day),
    then calls the `plot_shadow` function with the 'HOUR' mode to plot shadows 
    for the specified date, time, and location.
    """
    date = datetime(2025, 1, 1) + timedelta(days=day_year - 1)
    month, day = date.month, date.day
    plot_shadow(month, day, hour, minute, latitude_rad, longitude_rad, ax, points, 'HOUR')

def plot_evolution_location_in_which_sunrays_orthogonal(i, ax):
    """
    Plot the location where the Sun's rays are orthogonal (subsolar point) on the Earth's map
    for a specific hour `i` of the year.

    Parameters:
        i (int): Index representing the hour since Jan 1st, 00:00 UTC (0 to 8759)
        ax (matplotlib axis): The axis with a cartopy map projection to draw on
    """

    # Clear the previous frame
    ax.clear()

    # Add Earth map features
    land = cfeature.NaturalEarthFeature(
        'physical', 'land', '110m',
        edgecolor='face',
        facecolor='lightgreen'
    )
    ax.add_feature(land)
    ax.coastlines()
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.BORDERS)
    ax.set_global()

    # Add parallels and meridians
    ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')

    # Tropic of Cancer
    ax.plot(
        [-180, 180], [EARTH_ROTATION_AXIS_ANGLE_DEG, EARTH_ROTATION_AXIS_ANGLE_DEG],
        transform=ccrs.PlateCarree(), linestyle='--', linewidth=1, color='gray'
    )

    # Arctic Circle
    ax.plot(
        [-180, 180], [90 - EARTH_ROTATION_AXIS_ANGLE_DEG, 90 - EARTH_ROTATION_AXIS_ANGLE_DEG],
        transform=ccrs.PlateCarree(), linestyle='--', linewidth=1, color='gray'
    )

    # Tropic of Capricorn
    ax.plot(
        [-180, 180], [-EARTH_ROTATION_AXIS_ANGLE_DEG, -EARTH_ROTATION_AXIS_ANGLE_DEG],
        transform=ccrs.PlateCarree(), linestyle='--', linewidth=1, color='gray'
    )

    # Antarctic Circle
    ax.plot(
        [-180, 180], [-90 + EARTH_ROTATION_AXIS_ANGLE_DEG, -90 + EARTH_ROTATION_AXIS_ANGLE_DEG],
        transform=ccrs.PlateCarree(), linestyle='--', linewidth=1, color='gray'
    )

    # Determine the current day and hour
    day = i // 24
    hour = i % 24

    # Compute the subsolar point (latitude, longitude) for that specific hour
    (latitude_deg, longitude_deg) = location_in_which_sunrays_orthogonal(day + hour / 24)

    # Plot the subsolar point as a red dot
    ax.plot(longitude_deg, latitude_deg, 'ro', transform=ccrs.Geodetic())

    # Construct the date string for labeling
    date = datetime(2025, 1, 1) + timedelta(days=day)
    month, day = date.month, date.day

    # Create and set the plot title with date and subsolar point info
    date_label = MONTHS[month - 1] + ' ' + str(day) + ' UTC ' + str(hour).rjust(2, '0') + ':' + str(0).rjust(2, '0') + '\n'
    latitude_label = "Latitud: " + f"{latitude_deg:.2f}°".rjust(7)
    longitude_label = " Longitude: " + f"{longitude_deg:.2f}°".rjust(9)
    plt.title(date_label + latitude_label + longitude_label)


def generate_shadow_gifs_all_dates_all_latitudes():
    """
    Generate and save GIF animations of daily shadow patterns for multiple latitudes and selected dates.
    """
    figure, ax = plt.subplots()
    dates = [(2,3), (3,20), (5,4), (6,21), (8,6), (9,22), (11,6), (12,21)] # Solstices, equinoxes, and days in between
    for i in range(1,18): # 80°, 70°, ..., -70°, -80°
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

def generate_shadow_gif_given_date(latitude_deg, longitude_deg, month, day):
    """
    Generate and save a GIF animation showing the shadow pattern over a full day
    (24 hours) at a specific geographic location and date.

    Parameters:
    - latitude_deg: Latitude in degrees.
    - longitude_deg: Longitude in degrees.
    - month: Month of the year (1-12).
    - day: Day of the month.
    """
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

def generate_shadow_hour_gifs_all_latitudes():
    """
    Generate and save GIF animations of shadow patterns throughout the year
    at specific hours for multiple latitudes
    """
    figure, ax = plt.subplots()
    minute = 0
    for i in range(1,18): # 80°, 70°, ..., -70°, -80°
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

def generate_shadow_gif_given_hour(latitude_deg, longitude_deg, hour, minute):
    """
    Generate and save a GIF animation showing the shadow pattern over a full year
    at a specified hour and minute for a given geographic location.

    Parameters:
    - latitude_deg: Latitude in degrees.
    - longitude_deg: Longitude in degrees.
    - hour: Hour of the day (0-23) to visualize shadows.
    - minute: Minute of the hour (0-59) to visualize shadows.
    """

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

def plot_sunrise_and_sunset_times(latitude_rad):
    """
    Compute and plot daily sunrise and sunset times over one year for a given latitude.

    Parameters:
    - latitude_rad: Latitude in radians for which to calculate sunrise and sunset times.
    """
    sunrise=[]
    sunset=[]
    # For each day of the year 
    for day in range(365):

        # Define a function `f(t)` that represents the dot product between the
        # position vector of the location at fractional time `year_day + t` and 
        # the sun's vector projected onto the Earth's surface plane.

        f = lambda t:  np.dot(np.array(position(latitude_rad, 0, day + t)),
                               np.array((-math.cos(angle_sun_earth(day + t)), -math.sin(angle_sun_earth(day + t)), 0)))
        # Use the bisection method to find:
        #    * Sunrise time: root of f(t) between midnight and noon (0 to 0.5 fraction of day).
        #    * Sunset time: root of f(t) between noon and next midnight (0.5 to 1 fraction of day).
        # ... and add them to the list of sunrise and sunsets
        midnight = 0
        noon = 0.5 
        midnight_next_day = 1
        sunrise.append(bisection_method(midnight, noon, f, f(midnight), f(noon), 1/(24*60)))
        sunset.append(bisection_method(noon, midnight_next_day, f, f(noon), f(midnight_next_day), 1/(24*60)))
    
    # Create a list of datetime objects corresponding to each day in the year 2025.
    days = 365
    start_date = datetime(2025, 1, 1)
    dates = [start_date + timedelta(days=i) for i in range(days)]

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot sunrise and sunset times as lines against dates
    ax.plot(dates, sunset, label='Sunset', color='purple', linewidth=2)
    ax.plot(dates, sunrise, label='Sunrise', color='orange', linewidth=2)

    # Mark solstices and equinoxes
    event_dates = {
        "Spring Equinox": datetime(2025, 3, 20),
        "Summer Solstice": datetime(2025, 6, 21),
        "Autumn Equinox": datetime(2025, 9, 22),
        "Winter Solstice": datetime(2025, 12, 21)
    }

    for label, date in event_dates.items():
        ax.axvline(date, color='gray', linestyle='--', linewidth=1)
        ax.text(date + timedelta(days=-5), 0.4, label, rotation=90, va='bottom', ha='center', fontsize=9, color='gray')
    
    # Format the x-axis to show months, and the y-axis to show time of day in UTC.
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))

    yticks = np.linspace(0, 1, 25)  # 0 to 24h in 1-hour steps
    ytick_labels = [f"{int(hours):02d}:{int(0):02d}" for hours in range(25)]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels)

    # Titles and labels
    ax.set_xlabel('Month', fontsize=12)
    ax.set_ylabel('Time of Day (UCT)', fontsize=12)
    plt.suptitle('Sunrise and Sunset Times Over the Year\n\n'+'Latitude: ' + str(round(geom.radian2degree(latitude_rad),2)) + '°') 

    # Aesthetics: add grid, legend, and titles including the latitude.
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend()
    fig.tight_layout()

    plt.show()
 
def get_sunrise_and_sunset_time(latitude_rad, month, day):
    """
    Calculate and print the sunrise and sunset times (in UTC) for a given location and date.

    Parameters:
    - latitude_rad: Latitude of the location in radians.
    - month: Month of the year (1-12).
    - day: Day of the month.

    Notes:
    - The precision for the bisection method is set to 1 minute (1/(24*60)).
    """
    # Convert the given date to a fractional day of the year
    year_day = get_time_in_days(month, day, 0, 0)
    # Define a function `f(t)` that represents the dot product between the
    # position vector of the location at fractional time `year_day + t` and 
    # the sun's vector projected onto the Earth's surface plane.
    f = lambda t:  np.dot(np.array(position(latitude_rad, 0, year_day + t)),
                            np.array((-math.cos(angle_sun_earth(year_day + t)), -math.sin(angle_sun_earth(year_day + t)), 0)))
    
    # Use the bisection method to find:
    #    * Sunrise time: root of f(t) between midnight and noon (0 to 0.5 fraction of day).
    #    * Sunset time: root of f(t) between noon and next midnight (0.5 to 1 fraction of day).
    midnight = 0
    noon = 0.5 
    midnight_next_day = 1

    sunrise = bisection_method(midnight, noon, f, f(midnight), f(noon), 1/(24*60))
    sunset = bisection_method(noon, midnight_next_day, f, f(noon), f(midnight_next_day), 1/(24*60))

    # Convert fractional day times for sunrise and sunset into hours and minutes and 
    # print formatted sunrise and sunset times in UTC.
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

def generate_evolution_location_in_which_sunrays_orthogonal():
    """
    Generate and save an animated GIF showing how the location on Earth
    where the Sun's rays are orthogonal (subsolar point) changes throughout the year.

    The animation uses a Plate Carrée projection to display Earth's surface,
    and updates hourly over an entire year (365 days × 24 hours = 8760 frames).
    """

    # Create a matplotlib figure and axis with a Plate Carrée map projection,
    # which maps latitude and longitude directly to x and y coordinates.
    figure, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

    ani = animation.FuncAnimation(
        figure,
        func=plot_evolution_location_in_which_sunrays_orthogonal,
        frames=[i for i in range(365 * 24)],  # One frame per hour for 365 days
        fargs=(ax,),
        init_func=init_function
    )

    # Save the resulting animation as a GIF file
    ani.save('Orthogonal light.gif', writer='pillow', fps=1000)
          
def get_location_in_which_sunrays_orthogonal(month, day, hour, minute):
    """
    Compute and print the location (latitude and longitude) on Earth where the Sun's rays
    are orthogonal at a specific date and time in UTC.

    Parameters:
        month (int): Month of the year (1–12)
        day (int): Day of the month 
        hour (int): Hour in 24-hour format (0–23)
        minute (int): Minute (0–59)
    """

    # Convert the given date and time into a continuous time value `t` in days
    t = get_time_in_days(month, day, hour, minute)

    # Use the computed time to get the latitude and longitude where the sunrays are orthogonal
    (latitude_deg, longitude_deg) = location_in_which_sunrays_orthogonal(t)

    # Print the results
    print(MONTHS[month - 1] + ' ' + str(day) + ' UTC ' + str(hour).rjust(2, '0') + ':' + str(minute).rjust(2, '0'))
    print(f"Latitude: {latitude_deg:.2f}°, Longitude: {longitude_deg:.2f}°")




if __name__ == "__main__":
    mode = int(input(
        "Select a functionality by entering its corresponding number:\n"
        "1. Generate sun trajectory and shadow evolution during the day for solstices, equinoxes, and key in-between dates at informative latitudes.\n"
        "2. Generate sun trajectory and shadow evolution on a specific day at a given location.\n"
        "3. Generate sun and shadow evolution throughout the year at selected hours and informative latitudes.\n"
        "4. Generate sun and shadow evolution throughout the year at a specific time of day and location.\n"
        "5. Plot sunrise and sunset times over the year for a given latitude.\n"
        "6. Get sunrise and sunset times on a specific day at a given latitude.\n"
        "7. Generate the annual evolution of the location where sunrays are orthogonal to Earth’s surface (subsolar point).\n"
        "8. Get the location where sunrays are orthogonal to Earth’s surface for a specific date and time.\n"
        "Enter your selection (1–8): "
    ))

    if mode == 1:
        generate_shadow_gifs_all_dates_all_latitudes()
    elif mode == 2:
        latitude_deg = float(input("Enter the latitude (in decimal degrees, e.g., 40.4170): "))
        longitude_deg = float(input("Enter the longitude (in decimal degrees, e.g., -3.7034): "))
        month = int(input("Enter the month as a number (e.g., January = 1): "))
        day = int(input("Enter the day of the month: "))
        generate_shadow_gif_given_date(latitude_deg, longitude_deg, month, day)
    elif mode == 3:
        generate_shadow_hour_gifs_all_latitudes()
    elif mode == 4:
        latitude_deg = float(input("Enter the latitude (in decimal degrees, e.g., 40.4170): "))
        longitude_deg = float(input("Enter the longitude (in decimal degrees, e.g., -3.7034): "))
        hour = int(input("Enter the hour of the day (0–23): "))
        minute = int(input("Enter the minute (0–59): "))
        generate_shadow_gif_given_hour(latitude_deg, longitude_deg, hour, minute)
    elif mode == 5:
        latitude_deg = float(input("Enter the latitude (in decimal degrees, e.g., 40.4170): "))
        plot_sunrise_and_sunset_times(geom.degree2radian(latitude_deg))
    elif mode == 6:
        latitude_deg = float(input("Enter the latitude (in decimal degrees, e.g., 40.4170): "))
        month = int(input("Enter the month as a number (e.g., January = 1): "))
        day = int(input("Enter the day of the month: "))
        get_sunrise_and_sunset_time(geom.degree2radian(latitude_deg), month, day)
    elif mode == 7:
        generate_evolution_location_in_which_sunrays_orthogonal()
    elif mode == 8:
        month = int(input("Enter the month as a number (e.g., January = 1): "))
        day = int(input("Enter the day of the month: "))
        hour = int(input("Enter the hour of the day (0–23): "))
        minute = int(input("Enter the minute (0–59): "))
        get_location_in_which_sunrays_orthogonal(month, day, hour, minute)
