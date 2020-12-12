import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objs as go
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
from tqdm import tqdm
from plotly.offline import init_notebook_mode, iplot
import seaborn as sns
import plotly.figure_factory as ff


class Data_Explorer:
    
    def __init__(self, plotly=True):
        self.plotly=plotly
        

    #Beta value distribution (Illumina Analisys)
    def beta_distribution(self, dataset, labels):
        #select 4 samples for each label

        #TCGA-LUAD
        if not self.plotly:
            fig, axs = plt.subplots(4)
            fig.suptitle('TCGA-LUAD beta value distribution')
        else:
            fig = make_subplots(rows=1, cols=4)

        TCGA_LUAD_df = dataset.loc[labels=='TCGA-LUAD'].sample(4)
        i=0
        for index,row in TCGA_LUAD_df.iterrows():
            sample = row.array
            betas_in_sample={}
            for beta_value in sample:
                value="{:.2f}".format(float(beta_value))
                if float(value)<0.5:
                    if value in betas_in_sample:
                        betas_in_sample[value]+=1
                    else:
                        betas_in_sample[value]=1
            if not self.plotly:
                axs[i].bar(betas_in_sample.keys(),betas_in_sample.values())
            else:
                fig.add_trace(go.Bar(y=list(betas_in_sample.values())), row=1, col=i+1)
            i+=1
        if not self.plotly:
            plt.show();
        else:
            fig.update_layout(height=600, width=600, title_text="TCGA-LUAD beta value distribution")
            fig.show()

        #TCGA-LUSC
        if not self.plotly:
            fig, axs = plt.subplots(4)
            fig.suptitle('TCGA-LUSC beta value distribution')
        else:
            fig = make_subplots(rows=1, cols=4)

        TCGA_LUSC_df = dataset.loc[labels=='TCGA-LUSC'].sample(4)
        i=0
        for index,row in TCGA_LUSC_df.iterrows():
            sample = row.array
            betas_in_sample={}
            for beta_value in sample:
                value="{:.2f}".format(float(beta_value))
                if float(value)<0.5:
                    if value in betas_in_sample:
                        betas_in_sample[value]+=1
                    else:
                        betas_in_sample[value]=1
            if not self.plotly:
                axs[i].bar(betas_in_sample.keys(),betas_in_sample.values())
            else:
                fig.add_trace(go.Bar(y=list(betas_in_sample.values())), row=1, col=i+1)
            i+=1
        if not self.plotly:
            plt.show();
        else:
            fig.update_layout(height=600, width=600, title_text="TCGA-LUSC beta value distribution")
            fig.show()

        #False
        if not self.plotly:
            fig, axs = plt.subplots(4)
            fig.suptitle('Non tumoral beta value distribution')
        else:
            fig = make_subplots(rows=1, cols=4)

        Non_tumoral_df = dataset.loc[labels==False].sample(4)
        i=0
        for index,row in Non_tumoral_df.iterrows():
            sample = row.array
            betas_in_sample={}
            for beta_value in sample:
                value="{:.2f}".format(float(beta_value))
                if float(value)<0.5:
                    if value in betas_in_sample:
                        betas_in_sample[value]+=1
                    else:
                        betas_in_sample[value]=1
            if not self.plotly:
                axs[i].bar(betas_in_sample.keys(),betas_in_sample.values())
            else:
                fig.add_trace(go.Bar(y=list(betas_in_sample.values())), row=1, col=i+1)
            i+=1
        if not self.plotly:
            plt.show();
        else:
            fig.update_layout(height=600, width=600, title_text="Non tumoral beta value distribution")
            fig.show()


    #count per milion distribution (miRNA Analisys)
    def count_per_milion_distribution(self, dataset, labels):
        #TCGA-LUAD
        if not self.plotly:
            fig, axs = plt.subplots(4)
            fig.suptitle('TCGA-LUAD counts per milion value distribution')
        else:
            fig = make_subplots(rows=1, cols=4)

        TCGA_LUAD_df = dataset.loc[labels=='TCGA-LUAD'].sample(4)

        patients=[]
        for index, sample in TCGA_LUAD_df.iterrows():
            patients.append(sample.to_list())

        x_coordinates=[0]
        for i in range(0, len(patients[0])-1):
            x_coordinates.append(x_coordinates[i]+0.5)

        for i, patient in enumerate(patients):
            if not self.plotly:
                axs[i].bar(x_coordinates,patient)
            else:
                fig.add_trace(go.Bar(y=list(patient)), row=1, col=i+1)
            i+=1

        if not self.plotly:
            plt.show();
        else:
            fig.update_layout(height=600, width=600, title_text="TCGA-LUAD beta value distribution")
            fig.show()

        #TCGA-LUSC
        if not self.plotly:
            fig, axs = plt.subplots(4)
            fig.suptitle('TCGA-LUSC counts per milion value distribution')
        else:
            fig = make_subplots(rows=1, cols=4)

        TCGA_LUAD_df = dataset.loc[labels=='TCGA-LUSC'].sample(4)

        patients=[]
        for index, sample in TCGA_LUAD_df.iterrows():
            patients.append(sample.to_list())

        x_coordinates=[0]
        for i in range(0, len(patients[0])-1):
            x_coordinates.append(x_coordinates[i]+0.5)

        for i, patient in enumerate(patients):
            if not self.plotly:
                axs[i].bar(x_coordinates,patient)
            else:
                fig.add_trace(go.Bar(y=list(patient)), row=1, col=i+1)
            i+=1

        if not self.plotly:
            plt.show();
        else:
            fig.update_layout(height=600, width=600, title_text="TCGA-LUSC beta value distribution")
            fig.show()

        #False
        if not self.plotly:
            fig, axs = plt.subplots(4)
            fig.suptitle('Non tumoral counts per milion value distribution')
        else:
            fig = make_subplots(rows=1, cols=4)

        TCGA_LUAD_df = dataset.loc[labels==False].sample(4)

        patients=[]
        for index, sample in TCGA_LUAD_df.iterrows():
            patients.append(sample.to_list())

        x_coordinates=[0]
        for i in range(0, len(patients[0])-1):
            x_coordinates.append(x_coordinates[i]+0.5)

        for i, patient in enumerate(patients):
            if not self.plotly:
                axs[i].bar(x_coordinates,patient)
            else:
                fig.add_trace(go.Bar(y=list(patient)), row=1, col=i+1)
            i+=1

        if not self.plotly:
            plt.show();
        else:
            fig.update_layout(height=600, width=600, title_text="Non tumoral beta value distribution")
            fig.show()


    def labels_distribution(self, labels, title):
        distribution = [0, 0, 0]

        for label in labels:
            if label=='TCGA-LUAD':
                distribution[0] += 1

            if label=='TCGA-LUSC':
                distribution[1] += 1

            if label==False:
                distribution[2] += 1

        x_labels = ['TCGA-LUAD', 'TCGA-LUSC', 'Non Tumoral']
        fig = go.Figure(data=[go.Bar(x=x_labels, y=distribution, text=distribution, 
        marker_color=['rgb(26,35,126)', 'rgb(66,179,213)', 'rgb(220,237,200)'], textposition='auto',)])
        fig.update_layout(title_text=title, height=500, width=500)
        fig.show()
    
    def boxPlot(self, X, y, plotly=True):
        dataset = pd.concat([X, y], axis=1).T
        if plotly is True:
            
            fig = px.box(dataset, x="case_id", color=y)

            #fig.update_layout("Box plot",
            #                 boxmode='group',
            #                 xaxis=dict(title="case_id"))
            
            iplot(fig, filename="box-plot")
        else:
            ax = sns.boxplot(data=dataset, x="case_id", hue=y)
            ax.title("Box plot")
            plt.show()
        return
    
    def heatmap(self, X, y, plotly=True):
        correlation_matrix = X.corr()
        if plotly is True:
            trace = go.Heatmap(z=correlation_matrix.values.tolist(),
                              x=correlation_matrix.columns,
                              y=correlation_matrix.columns)
            data = [trace]
            layout=go.Layout(title="Correlation matrix")
            fig_heatmap = go.Figure(data=data, layout=layout)
            iplot(fig_heatmap)
        else:
            sns.color_palette("rocket")
            ax = sns.heatmap(correlation_matrix)
            ax.set_title("Correlation Matrix")
            plt.show()
        return
    
    def p_valueDistribution(self, dataset, for_tumor=False, plotly=False):
        if plotly is True:
            if for_tumor is True:
                group_labels = list(set(dataset.loc[:, 'tumor']))
                data = [dataset[dataset['tumor'] == l].loc[:, 'p_value'] for l in group_labels]
                group_labels = [str(i) for i in group_labels]
                fig = ff.create_distplot(data, group_labels)
                fig.show()
            else:
                group_labels = ['p_value']
                data = [dataset.loc[:, 'p_value'].values]
                fig = ff.create_distplot(data, group_labels)
                fig.show()
        else:
            if for_tumor is True:
                sns.distplot(data=dataset, x='p_value', hue='tumor')
                plt.show()
            else:
                sns.distplot(data=dataset)
                plt.show()
        return