'''
This piece of code is used for downloading altitude data and organizing altitude data into the final_summary
table. Final summary table is generated by concating all the pollutant tables together.

[Update: 01/26/2016] I added teh code for finding states using googleAPI
'''
#For appending to the final_summary dataframe. Final summary dataframe is the name of all
#the pollutant appended to one another.
#Loops through each row
for index,rows in final_summary.iterrows():
    a = rows['lat']
    b = rows['lon']
    #Checks whether there is a value in df_latitude that matches the longitude and latitude final_summary
    ss = df_altitude[(df_altitude.lat == a) & (df_altitude.lon == b)].values
    if (len(ss) != 0):
    	#If it is a number, return the value to the cell
        ss1 = ss[0][2]
        ss2 = ss[0][3]
        final_summary.set_value(index,'alt',ss1 )
        final_summary.set_value(index,'state',ss2)


#For doownloading the altitude data from google API
upper_lat = 50.0
lower_lat = 17.5
upper_lon = 310.0
lower_lon = 170.0
n_lat = 14
n_lon = 57

#declare empty dataframe
index = range(0,14*57)
df_altitude = pd.DataFrame(index = index, columns = ['lon','lat','alt','state'])
index=0

for lon in np.linspace(lower_lon,upper_lon,n_lon):
    for lat in np.linspace(lower_lat,upper_lat,n_lat):
        df_altitude.iloc[index,1] = lat
        if lon > 180:
            lon = lon-360 #return negative longitude, which is within domain
        url = 'https://maps.googleapis.com/maps/api/elevation/json?locations=' + str(lat)+','+str(lon)+'&key=AIzaSyCUaQB8Ilwyi8D0KGntFoa7_u3kYAjOAbg'
        df_altitude.iloc[index,0] = lon
        request = Request(url)
        response = urlopen(request)
        elevations = response.read()
        data = json.loads(elevations)
        el = []
        for result in data['results']:
            el.append(result[u'elevation'])
        #If there is no returned elevantion, then automatically set the elevantion to be zero
        if(len(el) != 0):
            df_altitude.iloc[index,2] = el[0]

        url2 = 'https://maps.googleapis.com/maps/api/geocode/json?latlng=' + str(lat)+','+str(lon)+'&key=AIzaSyCUaQB8Ilwyi8D0KGntFoa7_u3kYAjOAbg'
        request2 = urllib2.Request(url2)
        response2 = urllib2.urlopen(request2)
        geocode = response2.read()
        data_geocode = json.loads(geocode)
        if (data_geocode['results'] != null):
            length = len(data_geocode['results'])
            if (str(data_geocode['results'][length-1]['address_components'][0]['short_name']) == 'US'):
                state = str(data_geocode ['results'][length-2]['address_components'][0]['short_name'])
                df_altitude.iloc[index,3] = state
        index+=1


df_altitude = df_altitude[df_altitude.lon < 170]
df_altitude = df_altitude[df_altitude.lon < -65]

df_altitude = df_altitude[df_altitude.lat != 17.5]
df_altitude = df_altitude[df_altitude.lat != 20]
df_altitude = df_altitude[df_altitude.lat != 22.5]
