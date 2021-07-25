"""
Explore the SSD curves estimated in the antibiotic resistance project.
"""

import mizani
import pandas
import plotnine
import streamlit


MSC_FILEPATH = 'MSC.csv'
CURVEDATA_FILEPATH = 'curvedata.csv'
HC_FILEPATH = 'HC.csv'

CONCENTRATION_VALS = [2**x for x in range(-27, 6, 2)]
CONCENTRATION_LABELS = [str(round(x, 3)) if x > 0.001 else '' for x in CONCENTRATION_VALS]

TEXT_SLIDER_RANGE = (5, 15)
TEXT_SLIDER_STARTPOINT = 5
LINE_START = -27
LINE_FADE_FACTOR = 0.2
HC_LINETYPE = 'dotted'


# %% streamlit config

#streamlit.set_page_config(layout = 'wide')


# %% load data

@streamlit.cache
def load_data(*paths):

    return [pandas.read_csv(x) for x in paths]


d, curves, hc = load_data(
    MSC_FILEPATH,
    CURVEDATA_FILEPATH,
    HC_FILEPATH
)


# %% sidebar

antibiotic = streamlit.sidebar.selectbox(
    'antibiotic',
    pandas.unique(d['Antibiotic'])
)

d = d[d['Antibiotic'] == antibiotic]
curves = curves[curves['Antibiotic'] == antibiotic]
hc = hc[hc['Antibiotic'] == antibiotic]

casrn = d['CASRN'].iloc[0]


level = streamlit.sidebar.selectbox(
    'estimation unit',
    ['species', 'genus']
)

d = d[d['level'] == level]
hc = hc[hc['level'] == level]


show_names = streamlit.sidebar.checkbox(
    'show ' + level + ' names',
    value = True
)


if show_names:
    show_genus = streamlit.sidebar.checkbox('color by genus')
else:
    show_genus = False


if show_names:
    text_size = streamlit.sidebar.slider(
        'text size',
        min_value = TEXT_SLIDER_RANGE[0],
        max_value = TEXT_SLIDER_RANGE[1],
        value = TEXT_SLIDER_STARTPOINT,
        format = ''
    )


show_hc = []

for cutoff in pandas.unique(hc['cutoff']):
    if streamlit.sidebar.checkbox('show HC' + str(int(cutoff * 100))):
        show_hc.append(cutoff)

hc = hc[hc['cutoff'].isin(show_hc)]


# %% main column

streamlit.title(antibiotic + ' (CASRN: ' + casrn + ')')
streamlit.markdown('estimated by **' + level + '**')


fig = (
    plotnine.ggplot(
        d,
        plotnine.aes(x = 'MSC', y = 'PAF')
    )
    + plotnine.scale_x_continuous(
        breaks = CONCENTRATION_VALS,
        labels = CONCENTRATION_LABELS,
        trans = mizani.transforms.log2_trans
    )
    + plotnine.theme(
        axis_text_x = plotnine.element_text(rotation = 45, hjust = 1)
    )
    + plotnine.guides(color = False)
    + plotnine.labs(
        x = 'concentration (mg/L)',
        y = 'PAF ' + level
    )
    + plotnine.geom_segment(
        plotnine.aes(xend = 'HC', y = 'cutoff', yend = 'cutoff'),
        data = hc,
        x = LINE_START,
        linetype = HC_LINETYPE
    )
    + plotnine.geom_segment(
        plotnine.aes(x = 'HC', xend = 'HC', yend = 'cutoff'),
        data = hc,
        y = 0,
        linetype = HC_LINETYPE
    )
    + plotnine.geom_label(
        plotnine.aes(x = 'HC', y = 'cutoff', label = 'round(HC, 3)'),
        data = hc,
        ha = 'left'
    )
    + plotnine.geom_line(
            data = curves[curves['level'] != level],
            alpha = LINE_FADE_FACTOR
    )
    + plotnine.geom_line(
            data = curves[curves['level'] == level],
            color = 'black'
    )
)


if show_genus:
    names_aes = plotnine.aes(label = 'name', color = 'Genus')
else:
    names_aes = plotnine.aes(label = 'name')


if show_names:
    fig = (
        fig
        + plotnine.geom_text(
            names_aes,
            size = text_size,
            ha = 'right'
        )
    )


streamlit.pyplot(fig.draw())


# %% data (for debugging)

#streamlit.dataframe(d)
#streamlit.dataframe(curves)
#streamlit.dataframe(hc)
