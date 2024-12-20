# spatial-arrangement-POIs
# **Capturing the Spatial Arrangement of POIs in Crime Modeling**

We introduce two methods to capture the spatial arrangement of POIs: one is the Average Nearest Neighbor ratio (ANN_ratio), and another is the normalized Shannon Voronoi Diagram (n_SVDE)

**Abstract**: Many scholars have established that facilities represented by Points-of-Interests (POIs) may function as crime generators and attractors, influencing criminal activities. While existing measurements of POIs primarily rely on quantitative counts, this count-based approach overlooks the spatial arrangement of POIs within an area, which can also contribute to crime. This paper introduces two methods to capture the spatial arrangement characteristics of POIs. One is called the normalized Shannon Voronoi Diagram-based Entropy (n_SVDE). A Voronoi diagram is constructed based on the spatial distributions of POIs in an area, resulting in polygons, each corresponding to one POI. The area proportions of these polygons are then used to calculate Shannon Entropy. A low entropy value indicates a clustering pattern, while a high value reflects a dispersed distribution. The other is the average nearest neighbor distance ratio (ANN_ratio). It is a ratio of the average of the nearest distances of POIs in an area over the expected average. The effectiveness of these two methods is tested by using negative binominal models to explain street robberies in Cincinnati. Our findings show that the n_SVDE significantly explains street robbery, while the ANN_ratio shows no statistical significance. Specifically, a less clustered spatial distribution of POIs is positively associated with an increased likelihood of crime events, while a highly clustered distribution corresponds to a lower likelihood of crime. This study represents one of the pioneering implementations in explicitly examining the spatial configuration of POIs, contributing new insights into environmental criminology and providing valuable empirical evidence for enhancing place management and optimizing police patrols. 



**Keywords**: POI; Shannon entropy; Voronoi diagram; average nearest neighbor; crime


## **Paper**

If you find the paper useful, you may cite from here: 

Liu, L., Gu, X.*, Lan, M., Zhou, H., Chen, D., and Su, Z. (2025). Capturing the Spatial Arrangement of POIs in Crime Modeling. Computers, Environment and Urban Systems. 

@article{Gu@ComputationalCriminology,

  title={Capturing the Spatial Arrangement of POIs in Crime Modeling},
  
  author={Liu, Lin and Gu, Xin* and Lan, Minxuan and Zhou, Hanlin and Chen, Debao and Su, Zihan},
  
  journal={Computers, Environment and Urban Systems},
  
  volume={xx},
  
  number={xx},
  
  pages={xx},
  
  year={20xx},
  
  doi = {will update},
  
  publisher={Elsevier}
  
}

## **Data**
 Crime data: Cincinnati Police Department 
 
 Point-of-Interests (POIs): SafeGraph Inc., Cincinnati Metro Department
 
 Socioeconomic information: U.S. Census Bureau American Community survey 
