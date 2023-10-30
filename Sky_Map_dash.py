# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26  2023

@author: theperk
"""

# cd Documents\Astronomie
# python Dash_Sky_Map0d.py
#  http://127.0.0.1:8050/

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from skyfield.api import Star, load, wgs84, N, S, W, E
from skyfield.data import hipparcos, mpc, stellarium
from skyfield.constants import T0
from skyfield.projections import build_stereographic_projection
from skyfield import almanac

from datetime import datetime, date, timedelta
from pytz import timezone

from dash import Dash, html, dcc, Input, Output, callback, State
import dash_bootstrap_components as dbc


# https://viyaleta.medium.com/how-to-make-a-sky-map-in-python-a362bf722bb2
#https://rhodesmill.org/skyfield/examples.html#when-is-the-galactic-center-above-the-horizon

# Fichier des villes dont population > 1000 habitants (sinon trop de communes et ralentissement application)
fichier_communes = 'villes_pop_coord1000.csv'
df = pd.read_csv(fichier_communes, sep=';')

liste_communes = [ {'label' : str(df.loc[index, 'Nom_complet']) + ' ( ' + str(df.loc[index, 'nom_departement']) + ' )', 'value' : str(df.loc[index,'lat_degre']) + '/' + str(df.loc[index, 'long_degre'])} for index in df.index]

_COLUMN_NAMES = (
    'DSOID', 'RAdeg', 'DEdeg', 'Bmag', 'Vmag',
    'OType', 'MType', 'MajRarcmin', 'MinRarcmin', 'OAdegrees', 'RS', 'RSerror',
    'Plx', 'Plxerror', 'NRSdist', 'NRSdisterror', 'NGC', 'IC', 'M', 'C', 'B',
    'Sh2', 'VdB', 'RCW', 'LDN', 'LBN', 'Cr', 'Mel', 'PGC', 'UGC', 'Ced', 'Arp', 'VV', 'PK', 'PN',
    'SNR', 'ACO', 'HCG', 'ESO', 'VdBH', 'DWB', 'Tr', 'St', 'Ru', 'VdB-Ha',
)

nom_mois = ['Jan', 'Fév', 'Mars', 'Avr', 'Mai', 'Juin', 'Juil', 'Août', 'Sept', 'Oct', 'Nov' , 'Déc']

def dsos_load_dataframe(fobj):
    """
    Given an open file for `catalog.txt`, return a parsed dataframe.

    If your copy of ``catalog.txt`` has already been unzipped, pass the
    optional argument ``compression=None``.

    """
    try:
        from pandas import read_csv, set_option
    except ImportError:
        raise ImportError("NO PANDAS NO CANDO")

    fobj.seek(0)
    magic = fobj.read(2)
    compression = None #'gzip' if (magic == b'\x1f\x8b') else None
    fobj.seek(0)

    df = read_csv(
        fobj, sep='	', names=_COLUMN_NAMES, compression=compression,
        comment='#',
        usecols=['DSOID', 'RAdeg', 'DEdeg', 'Bmag', 'Vmag', 'M'],
        na_values=[''],
    )
    df.columns = (
        'dso_id', 'ra_degrees', 'dec_degrees', 'magnitudeB', 'magnitude',
        'messier_id',
    )
    df = df.assign(
        ra_hours = df['ra_degrees'] / 15.0,
        epoch_year = 2000.0,
    )
    df.loc[df['messier_id'] != 0.0, 'label'] = 'M' + df['messier_id'].astype(str)
    df = df[~df.messier_id.eq(0.0)]
    return df.set_index('dso_id')


zone = timezone('Europe/Paris')
now = zone.localize(datetime.now())
ts = load.timescale()
t = ts.from_datetime(now) #datetime(2023, 10, 26, 23, 0, 0)))
# 180 = South 0 = North
degrees = 0.0

# pour les heures de nuit noire
midnight = now.replace(hour=0, minute=0, second=0, microsecond=0)
next_midnight = midnight + timedelta(days=1)
t0 = ts.from_datetime(midnight)
t1 = ts.from_datetime(next_midnight)

# lieu pour démarrer
champi = wgs84.latlon(49.708217*N, 4.663512*E, elevation_m = 257).at(t)
#position = champi.from_altaz(alt_degrees = 90, az_degrees = degrees)

# Récupération des positions en provenance du JPL

eph = load('de421.bsp')
sun = eph['sun']
earth = eph['earth']

moon = eph['moon']

planets = ['Mercury', 'Mars Barycenter', 'Venus Barycenter', 'Jupiter Barycenter', 'Saturn Barycenter', 'Uranus Barycenter', 'Neptune Barycenter']
nom_planets = ['Mercure', 'Mars', 'Venus', 'Jupiter', 'Saturne', 'Uranus', 'Neptune']
eph_planets = []
for nom_p in planets:
    eph_planets.append(eph[nom_p])


# Récupération catalogue Hipparcos
with load.open(hipparcos.URL) as f:
    stardata = hipparcos.load_dataframe(f)
    
df_star_names = pd.read_csv('Hip_stars_names.csv', sep=';')

stardata = pd.merge(stardata, df_star_names, how = 'left', left_index = True, right_on = 'hip_number')
stardata = stardata.set_index('hip_number')
stardata['nom'] = stardata['nom'].fillna('')

# Ouverture fichier des DSOs

with open('catalog.txt') as f:
    dsodata = dsos_load_dataframe(f)


# Récupération des contours des constellations

file_cons = 'constellationship.fab'
with open(file_cons, 'rb') as f:
    consdata = stellarium.parse_constellations(f)
    
    
def nuit_noire(jour_midnight, location):
    # pour avoir l'éphéméride du jour des heures de nuit civil, nautique et astronomique
    
    zone = timezone('Europe/Paris')         
    next_midnight = jour_midnight + timedelta(days=1)

    ts = load.timescale()
    t0 = ts.from_datetime(jour_midnight)
    t1 = ts.from_datetime(next_midnight)
    eph = load('de421.bsp')
    lati, longi = map(float,location.split('/'))
    lat = N if lati > 0 else S
    lon = E if longi > 0 else W
    localisation = wgs84.latlon(lati * lat, 
                                longi * lon)     
    
    f = almanac.dark_twilight_day(eph, localisation)
    times, events = almanac.find_discrete(t0, t1, f)

    previous_e = f(t0).item()
    dico = {}
    for t, e in zip(times, events):
        tstr = str(t.astimezone(zone))[:16]
        if previous_e < e:
            #print(tstr, ' ', almanac.TWILIGHTS[e], 'starts')
            dico[almanac.TWILIGHTS[e]+ ' starts'] = int(''.join(tstr[-5:].split(':')))
        else:
            #print(tstr, ' ', almanac.TWILIGHTS[previous_e], 'ends')
            dico[almanac.TWILIGHTS[previous_e]+ ' ends'] = int(''.join(tstr[-5:].split(':')))
        previous_e = e  
    
    return dico


def couleur_ciel(heure, ephem):
    heure = int(''.join(heure.split(':')))
    
    color = "#111111"
    colorext = '#222222'
    ambiance = 'nuit astronomique'
    if 'Astronomical twilight starts' in ephem:
        if heure < ephem['Astronomical twilight starts']:
            return ('#111111', '#222222', 'nuit astronomique')
    if 'Nautical twilight starts' in ephem:
        if heure < ephem['Nautical twilight starts']:
            return ('#1b2134', '#242636', 'aube astronomique/nuit nautique')
    if 'Civil twilight starts' in ephem:
        if heure < ephem['Civil twilight starts']:
            return ('#1f2f55', '#2f3f4f', 'aube nautique/nuit civile')
    if 'Day starts' in ephem:
        if heure < ephem['Day starts']:
            return ('#1f2f99', '#1f4f99', 'aube civile')
    if 'Day ends' in ephem:
        if heure < ephem['Day ends']:
            return ('#f0dabd', '#b4af9f', 'jour') 
    if 'Civil twilight ends' in ephem:
        if heure < ephem['Civil twilight ends']:
            return ('#1f2f99', '#1f4f99','crépuscule civil')
    if 'Nautical twilight ends' in ephem:
        if heure < ephem['Nautical twilight ends']:
            return ('#1f2f55', '#2f3f4f', 'crépuscule nautique / nuit civile')    
    if 'Astronomical twilight ends' in ephem:
        if heure < ephem['Astronomical twilight ends']:
            return ('#1b2134', '#242636', 'crépuscule astronomique / nuit nautique')
    return (color, colorext, ambiance)
        
   
    return color
    

def generate_constellation_lines(data, polygon = False):
    edges = [edge for name, edges in data for edge in edges]
    noms = [name for name, edges in data for edge in edges]
    
    edges_star1 = [star1 for star1, star2 in edges]
    edges_star2 = [star2 for star1, star2 in edges]
    xy1 = stardata[['x', 'y']].loc[edges_star1].values
    xy2 = stardata[['x', 'y']].loc[edges_star2].values  

    if polygon:
        return [xy1]
    else:

        # The constellation lines will each begin at the x,y of one star and end
        # at the x,y of another.  We have to "rollaxis" the resulting coordinate
        # array into the shape that matplotlib expects.

        return np.rollaxis(np.array([xy1, xy2]), 1), noms


def trace_fig(t, localisation, choix_options, mag_lim, dso_lim, color_sky, color_externe):
    
    # Récupère position du lieu d'observation
    position = localisation.from_altaz(alt_degrees = 90, az_degrees = degrees)
    
    # Projection
    projection = build_stereographic_projection(position)
    field_of_view_degrees = 180.0
    limiting_magnitude = mag_lim #6.0
    dso_limit_magnitude = dso_lim #8.0

    # Calcul coordonnées par projection, pour les étoiles et les DSOs, la Lune, le Soleil

    star_positions = earth.at(t).observe(Star.from_dataframe(stardata))
    stardata['x'], stardata['y'] = projection(star_positions)

    dso_positions = earth.at(t).observe(Star.from_dataframe(dsodata))
    dsodata['x'], dsodata['y'] = projection(dso_positions)
    
    moon_positions = earth.at(t).observe(moon)
    x_moon, y_moon = projection(moon_positions)
    
    sun_positions = earth.at(t).observe(sun)
    x_sun, y_sun = projection(sun_positions)
    
    planet_positions = [earth.at(t).observe(planet) for planet in eph_planets]
    x_plan = []
    y_plan = []
    nom_plan = []
    for plan in planet_positions:
        x_p, y_p = projection(plan)
        x_plan.append(x_p)
        y_plan.append(y_p)
     
    
    # Conserve (par masque) uniquement les étoiles dont magnitude inférieure à celle choisie

    bright_stars = (stardata.magnitude <= limiting_magnitude)
    magnitude = stardata['magnitude'][bright_stars]
    star_size = (0.5 + limiting_magnitude - magnitude) ** 2.0
    star_size = magnitude.apply(lambda x : max(1, 12 - x)) 
    criteria = [magnitude.between(-30, 1),magnitude.between(1, 3), magnitude.between(3, 5), magnitude.between(5, 25)]
    values = [15, 10, 6, 1]

    star_size = np.select(criteria, values, 0)
    
    max_star_size = 15
    star_size = max_star_size * 10 ** (magnitude / -3.5)

    bright_dsos = (dsodata.magnitude <= dso_limit_magnitude)
    dso_magnitude = dsodata['magnitude'][bright_dsos]
    dso_size = (0.9 + dso_limit_magnitude - dso_magnitude) ** 2.0
    criteria = [dso_magnitude.between(1, 4), dso_magnitude.between(4, 6), dso_magnitude.between(6, 8), dso_magnitude.between(8, 10), dso_magnitude.between(10, 12), dso_magnitude.between(12, 30)]
    values = [12, 8, 6, 4, 2, 1]
    dso_size = np.select(criteria, values, 0)
    #dso_size = 20 - 1.5*dso_magnitude
    
    #max_dso_size = 30
    #dso_size = max_dso_size * 10 ** (dso_magnitude / -3.5)


    #fig = go.Figure(layout = layout)
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    #fig.update_layout(layout=layout)
    #fig.data = []
    fig.update_layout(autosize=True,
                        width = 750,
                        height = 750,
                        xaxis=go.layout.XAxis(linecolor="black", linewidth=1, mirror=True),
                        yaxis=go.layout.YAxis(linecolor="black", linewidth=1, mirror=True),
                        plot_bgcolor = color_externe,                     
                        showlegend = False,
                        margin = dict(l = 20, r = 20, t = 20, b = 20),
                        paper_bgcolor = "LightSteelBlue",
                      title = dict({'text' : 'N', 'x' : 0.5, 'xanchor': 'center', 'yanchor' : 'bottom'})
                     )

    # Tracé des lignes d'altitudes     
    if 'azimut' in choix_options :
        for alti in range(10,90,10):
            h0 = projection(localisation.from_altaz(alt_degrees = alti, az_degrees = 0.0))
            horizon_x = [h0[0]]
            horizon_y = [h0[1]]
            for i in range(1, 73):
                delta = 5.0
                current = i * delta
                h1 = projection(localisation.from_altaz(alt_degrees = alti, az_degrees = current))
                horizon_x.append(h1[0])
                horizon_y.append(h1[1])       

            fig.add_trace(go.Scatter(
                x = horizon_x,
                y = horizon_y,
                mode = "lines",
                line = go.scatter.Line(color = "grey", dash = 'dash'),
                showlegend = False))

        # Tracé des lignes d'azimuths
        for az in range(0,360,30):
            h0 = projection(localisation.from_altaz(alt_degrees = 0, az_degrees = az))
            horizon_x = [h0[0]]
            horizon_y = [h0[1]]
            for i in range(1, 17):
                delta = 5.0
                current = i * delta
                h1 = projection(localisation.from_altaz(alt_degrees = current, az_degrees = az))
                horizon_x.append(h1[0])
                horizon_y.append(h1[1])       

            fig.add_trace(go.Scatter(
                x = horizon_x,
                y = horizon_y,
                mode = "lines",
                line = go.scatter.Line(color = "grey", dash = 'dash'),
                showlegend = False))
        
    # Tracé de l'horizon  
    fig.add_shape(type = "circle", xref = "x", yref = "y", fillcolor = color_sky,  opacity = 0.3,  x0 = -1, y0 = -1, x1 = 1, y1 = 1, line_color = "LightSeaGreen")
    
    # Tracé des constellations
    if 'constel' in choix_options:
        retour, noms = generate_constellation_lines(consdata)        
        
        for nombre, nom in zip(retour, noms):
            hx = [nombre[0][0], nombre[1][0]]
            hy = [nombre[0][1], nombre[1][1]]
            fig.add_trace(go.Scatter(
                    x = hx,
                    y = hy,
                    text = [nom],
                    textfont_color = 'white',
                    mode = "lines",                  
                    
                    line = go.scatter.Line(color = "lightgrey", width = 1, dash = 'dash'),
                    hovertemplate='%{text}',
                    showlegend = False))
  
    
    # ajout de la Lune   
    if 'nom_planet' in choix_options:
        mode = 'markers+text'
    else:
        mode = 'markers'
    fig.add_trace(go.Scatter(x = [x_moon], y = [y_moon],
                             mode = mode,
                             text =  ['Lune'],
                             textfont_color = 'white',
                             textposition = 'bottom center',
                             textfont_size = 15,
                             hovertemplate='%{text}',
                             name ='', # pour supprimer le nom de la trace dans l'infobulle
                             marker = dict(size = 15, color = 'white')))
    # ajout du soleil   
    if 'nom_planet' in choix_options:
        mode = 'markers+text'
    else:
        mode = 'markers'
    fig.add_trace(go.Scatter(x = [x_sun], y = [y_sun],
                             mode = mode,
                             text =  ['Soleil'],
                             textfont_color = 'white',
                             textposition = 'bottom center',
                             textfont_size = 15,
                             hovertemplate='%{text}',
                             name ='', # pour supprimer le nom de la trace dans l'infobulle
                             marker = dict(size = 15, color = 'white')))
      
    # ajout des planètes   
    fig.add_trace(go.Scatter(x = x_plan, y = y_plan,
                             mode = mode,
                             text =  nom_planets,
                             textfont_color = 'white',
                             textfont_size = 11,
                             textposition = 'bottom center',
                             hovertemplate='%{text}',
                             name ='', # pour supprimer le nom de la trace dans l'infobulle
                             marker = dict(size = 10, color = 'white')))
    
    
    # Tracé des DSOs, en rouge 
    if 'nom_dso' in choix_options:
        mode = 'markers+text'
    else:
        mode = 'markers'
    if dso_magnitude.shape[0] >= 1 :
        fig.add_trace(go.Scatter(x = list(dsodata['x'][bright_dsos].values), y = list(dsodata['y'][bright_dsos].values),
                                 mode = mode,
                             text = dsodata['label'][bright_dsos] ,
                             meta = '<br>(magnitude ' + dsodata['magnitude'][bright_dsos].apply(lambda x :str(x)) + ')',
                             textposition = 'bottom center',
                             textfont_color = 'red',
                             hovertemplate = '%{text}%{meta}',                             
                             name = '',
                             marker = dict(size = dso_size, color = 'red'))) #dso_size
        
    # Tracé des étoiles
    if 'nom_etoile' in choix_options:
        mode = 'markers+text'        
    else:
        mode = 'markers'       
    
    fig.add_trace(go.Scatter( x = stardata['x'][bright_stars], y = stardata['y'][bright_stars],                            
                             mode = mode,
                             text =  stardata['nom'][bright_stars],
                             textposition = 'bottom center',
                             textfont_color = 'white',
                             meta = '<br>(magnitude ' + stardata['magnitude'][bright_stars].apply(lambda x: str(x)) + ')',
                             hovertemplate='%{text}%{meta}',
                             name ='', # pour supprimer le nom de la trace dans l'infobulle
                             marker = dict(size = star_size, color = 'white')))
      
   
        
    # Tracé final en tenant compte des limites de vue
    angle = np.pi - field_of_view_degrees / 360.0 * np.pi
    limit = np.sin(angle) / (1.0 - np.cos(angle))
    
    # Ajout des repères cardinaux
    
    fig.update_xaxes(range = (-limit, limit),
                     #showticklabels = False,
                     ticktext=["", "S", ""],
                     tickvals=["-1", "0", "1"],
                     tickfont=dict(color = 'black', size = 14),
                     showgrid = False,
                     zeroline = False,
                    )
    fig.update_yaxes(range = (-limit, limit),
                     #showticklabels = False,
                     showgrid = False,
                     ticktext=["", "E", ""],
                     tickvals=["-1", "0", "1"],
                     tickfont=dict(color = 'black', size = 14),
                     zeroline = False,
                     secondary_y = False
                    )
    fig.update_yaxes(range = (-limit, limit),
                     #showticklabels = False,
                     showgrid = False,
                     ticktext=["", "O", ""],
                     tickvals=["-1", "0", "1"],
                     tickfont=dict(color = 'black', size = 14),
                     zeroline = False,
                     secondary_y = True
                    )
    fig.add_trace(go.Scatter(x=[-1,0,1], y=[-1,0,1],
                             name="yaxis2 data",
                             mode = 'text',
                             text = ['','','']),
                             secondary_y=True,
                 )
    
    
    return fig


def trace_fig_bar_data(df_data, mag_lim, titre):
        
    if mag_lim > 1:
         # récupère étoiles dont magnitude inférieure à celle souhaitée et visibles au lieu d'observation (donc dans un rayon de 1 sur la carte)
        bright_stars = (df_data.magnitude <= mag_lim) & (df_data['x']**2 + df_data['y']**2 <= 1)
        df_data2 = df_data[bright_stars].copy()
        bins = list(range(0, int(mag_lim) + 1, 1))
        group_names = [str(n-1) + ' < Mag <= ' + str(n) + ' ' for n in range(1, int(mag_lim) + 1)]
        df_data2['mag_range'] = pd.cut(df_data2['magnitude'], bins, labels = group_names)
        df2 = df_data2.groupby('mag_range').count().reset_index()
        df2 = df2.rename(columns = {"magnitude" : "compte_mag"})

        fig3 = px.bar(df2, y = 'mag_range', x = 'compte_mag', text = 'compte_mag', orientation = 'h', color_discrete_sequence =['navy']*len(df2))
        
    else:
        fig3 = px.bar(pd.DataFrame({'x' : [0], 'y' : [0]}), y = "y", x = "x", orientation = 'h')
        fig3.update_yaxes(visible = False, showgrid= False)
        fig3.update_xaxes(title = '', visible = False)
        
    fig3.update_yaxes(title = '', visible = True)
    fig3.update_xaxes(title = '', visible = True)
    fig3.update_layout(title = dict({'text' : titre + ' ( ' + str(df2['compte_mag'].sum()) + ' )', 'x' : 0.5, 'xanchor': 'center', 'yanchor' : 'bottom'}),
                       margin = dict(l = 1, r = 1, t = 30, b = 10),
                       autosize = True,
                       paper_bgcolor = "LightSteelBlue")
    
    return fig3


fig_sky = trace_fig(t, champi, ['constel'], 5.0, 8.0, "#000023", "#000023")
fig_star = trace_fig_bar_data(stardata, 5.0, "Nombre d'étoiles\n ")
fig_dso = trace_fig_bar_data(dsodata, 8.0, "Nombre de DSO")

# liste des heures sélectionnables
heures = []
for h in range(24):
    for m in range(0,60,15):
        heures.append(str(h).zfill(2) + ':' + str(m).zfill(2))        
        
        
# DEBUT APPLICATION DASH -----------------

app = Dash(external_stylesheets = [dbc.themes.CYBORG]) #BOOTSTRAP])
server = app.server

app.layout = dbc.Container([    
    
    dbc.Row([
        dbc.Col(
                           
            dcc.Markdown('''            
               
                 Choix du lieu :
                
                 '''), style={'textAlign': 'center'},
                     width = 3),
        dbc.Col(
            dcc.Markdown('''                   
                 
              
                 
              Choix de la date :

                  '''),
            style = {'textAlign': 'center'},
            width = 2),
        
        dbc.Col(
            dcc.Markdown('''                   
                 
              
                 
              Choix de l'heure :

                  '''),
            style = {'textAlign': 'center'},
            width = 2),
        
        dbc.Col(
            dcc.Markdown('''                   
                 
              
                 
              Magnitude limite étoiles :

                  '''),
            width = 2),
        
        dbc.Col(
            dcc.Markdown('''                   
                 
              
                 
              Magnitude limite DSO :

                  '''),
            width = 2)
                ]),
                
    dbc.Row([
        dbc.Col(html.Div([
            # choix du lieu            
            dcc.Dropdown(id = "drop-location",
                         options = liste_communes,
                         value = liste_communes[531]['value'] #Charleville 
                        ),
            
                    ]),
                width = 3),
        
        dbc.Col(
            # choix de la date ;-))
            dcc.DatePickerSingle(id ='date-picker',
                                 month_format = 'MMMM YYYY',
                                 placeholder = 'DD MM YYYY',
                                 display_format ='DD MMM YYYY',
                                 min_date_allowed=date(1990, 1, 1),
                                 max_date_allowed=date(2030, 12, 12),
                                 initial_visible_month=date(2024, 10, 26),
                                 date=date.today(), #(2023, 8, 25)
                                 style={'width': '150px', 'margin-left': 30}
                                ),
            width = 2),
        
        dbc.Col(
           # choix de l'heure
                dcc.Dropdown(id = "heure-picker",
                             options = [{'label':  heure, 'value': heure} for heure in heures],
                             value = '23:00',
                             ),
                width = 2),
            dbc.Col(
                # choix mag lim étoiles
                dcc.Dropdown(id = "items_etoiles",
                               options = [{"label": f" Magnitude {nombre}", "value": nombre} for nombre in range(2, 11)],
                               value = 6
                              ),
                width = 2),
   
            dbc.Col(
                # choix mag lim DSO
                dcc.Dropdown(id = "items_dso",
                               options = [{"label" : " Pas de DSO", "value" : 0}] + [{"label": f" Magnitude {nombre}", "value": nombre} for nombre in range(4, 12)],
                               value = 8
                              ),
                width = 2)
    ]),
    dbc.Row(
            dbc.Col(
                html.H1(""),
                width = 12,
                style = {"height": "30px"},
            )
        ),
    dbc.Row([dbc.Col(
        # choix des options 
            dcc.Checklist(id = "choix-options",
                          options=[
                              {"label": " constellations", "value": "constel"},
                              {"label": " grille azimutale", "value": "azimut"},
                              {"label": " noms des planètes / Lune / Soleil", "value": "nom_planet"},
                              {"label": " noms des étoiles", "value": "nom_etoile"},
                              {"label": " noms des DSO", "value": "nom_dso"}                              
                                  ],
                          labelStyle = {"display" : "inline", 'margin-left' : '70px'},
                          value = ['constel', 'nom_planet'],
                         ), width = 12) 
    ]),
    dbc.Row(
            dbc.Col(
                html.H1(""),
                width = 12,
                style = {"height": "30px"},
            )
        ),
      dbc.Row([
            
        #dbc.Col(html.H4(children = 'Sky Map', id = 'titre_carte'), width = 8),
        dbc.Col(dcc.Textarea(value = 'Carte du ciel', id = 'titre_carte', style={'width': '100%', 'height': 58, 'background-color' : 'LightSteelBlue', 'color': 'navy'}), width = 7),
        dbc.Col(width = 1),
        dbc.Col(dcc.Textarea(value = "Nombre d'étoiles et de DSO (Objets du ciel profond)\nvisibles (par magnitude)", id = 'titre_bar_graph',style={'width': '100%', 'height': 58, 'background-color' : 'LightSteelBlue', 'color': 'navy'}), width = 4)
      ]),
    dbc.Row([
        dbc.Col(html.Div([
            dcc.Graph(id ='graph_sky',
                      figure = fig_sky,
                      className="flex-grow-1")],
            className="h-100 d-flex flex-column"),
            width = 8),  
        
        dbc.Col(html.Div([                    
            dcc.Graph(id ='graph_bar_stars',
                      figure = fig_star,
                      className="flex-grow-1"),
            dcc.Graph(id ='graph_bar_dsos',
                      figure = fig_dso,
                      className="flex-grow-1")],
            className="h-100 d-flex flex-column"),
                width = True),
        ])
         
    ])
                 
               
@callback(
    [Output('graph_sky', 'figure'),
     Output('titre_carte', 'value'),
     Output('graph_bar_stars', 'figure'),
     Output('graph_bar_dsos', 'figure')],
    [Input('drop-location', 'value'),
     Input('date-picker', 'date'),
     Input('heure-picker', 'value'),
     Input('choix-options', 'value'),     
     Input('items_etoiles', 'value'),
     Input('items_dso', 'value')
    ],
    [State("drop-location","options")]
)

def update_graph( location, date1, heure1, choix_options, maglim, dsolim, opt):
    
    lieu = [x['label'] for x in opt if x['value'] == location]
    retour = 'Carte du Ciel vue de ' +  str(lieu[0]) + '\n le ' + date1.split('-')[2] + ' ' + nom_mois[int(date1.split('-')[1])-1] + ' ' + date1.split('-')[0] + ' à ' + heure1 
    
    zone = timezone('Europe/Paris')
    ts = load.timescale()    
    t = ts.from_datetime(zone.localize(datetime(int(date1.split('-')[0]), int(date1.split('-')[1]), int(date1.split('-')[2]), int(heure1.split(':')[0]), int(heure1.split(':')[1]), 0)))#(2023, 10, 26, 23, 0, 0)))
    lati, longi = map(float,location.split('/'))
    lat = N if lati > 0 else S
    lon = E if longi > 0 else W
    localisation = wgs84.latlon(lati * lat, 
                                longi * lon,
                                elevation_m = 100).at(t) 
    jour_minuit = zone.localize(datetime(int(date1.split('-')[0]), int(date1.split('-')[1]), int(date1.split('-')[2]), 0, 0, 0))
    
    # récupère éphéméride du jour
    dico_ephem = nuit_noire(jour_minuit, location)
    
    # pour savoir quelle couleur de fond de ciel attribuer
    colorsky, colorexterne, ambiance = couleur_ciel(heure1, dico_ephem)
    #print(colorsky, colorexterne, ambiance)
    
    fig_sky = trace_fig(t, localisation, choix_options, maglim, dsolim, colorsky, colorexterne)
    fig_stars = trace_fig_bar_data(stardata, maglim, "Nombre d'étoiles\n ")
    fig_dsos = trace_fig_bar_data(dsodata, dsolim, "Nombre de DSO")
    
    return [fig_sky, retour + ' (' + ambiance + ')', fig_stars, fig_dsos]

    
if __name__ == '__main__':
    app.run(debug = True)