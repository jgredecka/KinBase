from flask import Flask, render_template, url_for, redirect
import pandas as pd

from flask_wtf import FlaskForm
from wtforms import StringField, SelectField
from wtforms.validators import Required

from sqlalchemy import Column, Integer, String, create_engine, ForeignKey, Table
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker, relationship, backref

app = Flask(__name__)
app.config['SECRET_KEY'] = '19c6553e-99ed-4ee8-84be-fed6c78d21bc'

# SQLAlchemy initial configuration:
# Define a type of connection to the database.
engine = create_engine('sqlite:///static/database/Howard_database_V4.db', echo=True)
# declarative_class creates a class with the engine. It is then used to create subclasses.
Base = declarative_base(engine)
metadata = Base.metadata
# Create a session so the database can be queried.
db_session = scoped_session(sessionmaker(autoflush=False, bind=engine))
global_session = db_session()

# Remove session at the end of request or when the app shuts down.
@app.teardown_appcontext
def shutdown_session(exception=None):
    db_session.remove()

# Define a class for the Kinase table which is mapped directly to the database. The class attributes correspond to columns of the table.
class Kinases(Base):
    __tablename__ = 'Kinase'
    
    Kinase_full_name = Column(String)
    Kinase_uniprot_code = Column(String)
    Gene_symbol = Column(String, primary_key=True)
    Kinase_family = Column(String)
    Cell_location = Column(String)
    inhibitors = relationship("Inhibitors", secondary="kin_inhi_relat", backref="Kinase")
    
    def __init__(self, Kinase_full_name, Kinase_uniprot_code, Gene_symbol, Kinase_family, Cell_location):

        self.Kinase_full_name = Kinase_full_name
        self.Kinase_uniprot_code = Kinase_uniprot_code
        self.Gene_symbol = Gene_symbol
        self.Kinase_family = Kinase_family
        self.Cell_location = Cell_location
    
# Define a class for the Inhibitor table in the database.
class Inhibitors(Base):
    __tablename__ = 'Inhibitor'
    
    Inhibitor = Column(String, primary_key=True)
    Chemical_Structure = Column(String)
    Molecular_Weight = Column(Integer)
    Chemical_image = Column(String)
    kinases = relationship("Kinases", secondary="kin_inhi_relat", backref="Inhibitor")
    
    def __init__(self, Inhibitor, Chemical_structure, Molecular_weight, Chemical_image):

        self.Inhibitor = Inhibitor
        self.Chemical_Structure = Chemical_structure
        self.Molecular_Weight = Molecular_weight
        self.Chemical_image = Chemical_image
    
# Define the association table for the many-to-many relationships.
kin_inhi_relat = Table('kin_inhi_relat', Base.metadata,
    Column('Inhibitor', String, ForeignKey(Inhibitors.Inhibitor), primary_key=True),
    Column('Kinases', String, ForeignKey(Kinases.Gene_symbol), primary_key=True)
    )

# Read the phosphosite data CSV file, convert it to a pandas dataframe and store as a global variable.
phospho_data = "static/files/final_phospho.csv"
read_df=pd.read_csv(phospho_data)

# Define a class using the imported FlaskForm. The variable stores a search form with relevant validators.
class QueryForm(FlaskForm):
    choices = [("Substrate", "Substrate Gene"), ("Kinase", "Kinase Gene"), ("Sequence (-)", "Sequence (-)"), ("Sequence (+)", "Sequence (+)"), ("Locus", "Substrate Gene Locus")]
    select = SelectField(choices=choices)
    user_query = StringField(validators=[Required()])
    
class KinInhForm(FlaskForm):
    choices = [("kinase", "Kinase"), ("inhibitor", "Inhibitor(s)")]
    select_mol = SelectField(choices=choices)
    search_term = StringField(validators=[Required()])

@app.route("/")
def index():
    return render_template('index.html', title='Home')

@app.route("/contact-us")
def contact():
    return render_template('contact.html', title='Contact Us')

# If the form is validated here, the user is redirected to the phosphosite results page.
@app.route("/browse", methods=['GET', 'POST'])
def prebrowse():
    form = QueryForm()
    user_query = None
    select = None
    if form.validate_on_submit():
        user_query = form.user_query.data
        select = form.select.data
        return redirect(url_for('browse', selection=select, query=user_query))
    return render_template('browse_phospho.html', title='Browse Phosphosites', form=form, user_query=user_query)

# Dataframe is manipulated according to user query:
# Rows whose column value equals the query are selected.
# Duplicates returned as a result of the function are dropped and hyperlinks are added in real time.
# The user is only redirected to the error message page if there are no results.
@app.route("/browse/<selection>/<query>")
def browse(selection, query):
    try:
        user_choice = selection 
        query=query.upper()
        df=read_df.drop_duplicates()
        df.set_index(user_choice, inplace=True, drop=False)
        df=df.loc[[query]]
        df['Kinase']=df['Kinase'].apply(lambda x: '<a href=/search/kinases/'+str(x)+'>'+str(x)+'</a>')
        number=len(df)
        return render_template("browse.html", title='Browse Phosphosites: ' + query, number=number, query=query, data=df.to_html(index=False, escape=False, classes = 'stripe hover" id = "myTable'))
    except KeyError:
        return redirect(url_for('user_message', query=query))

# If the kinase search form is validated here, the user is redirected to the kinase results page.
@app.route("/search", methods=['GET', 'POST'])
def search():
    #form = SearchForm()
    #search_term = None
    form = KinInhForm()
    select_mol = None
    search_term = None
    if form.validate_on_submit():
        select_mol = form.select_mol.data
        search_term = form.search_term.data
        search_term = search_term.upper()
        if select_mol == "kinase":
            return redirect(url_for('kinase_results', query=search_term))
        elif select_mol == "inhibitor":
            return redirect(url_for('inhibitor_list', query=search_term))
    return render_template("search.html", title='Search Kinases', form=form, search_term=search_term, select_mol=select_mol)

# A new session is loaded and kinase matching the query stored in the 'kinase' variable.
# Instance attributed are used to extract background info on the kinase.
@app.route("/search/kinases/<query>")
def kinase_results(query):
    kinase = global_session.query(Kinases).get(query)
    if kinase is None:
         return redirect(url_for('user_message', query=query))
    name = kinase.Kinase_full_name
    pr_gene = kinase.Gene_symbol
    families = kinase.Kinase_family
    families = families.rstrip(".")
    loc = kinase.Cell_location
    return render_template("kinase_results.html", title='Search Kinases: ' + query, query=query, name=name, gene=pr_gene, family=families, location=loc)

# New session is loaded and inhibitor(s) corresponding to the query are stored in 'inhibitor_list'.
# If inhibitors were found, the user is taken to the inhibitor search results page.
@app.route("/search/inhibitors/<query>")
def inhibitor_list(query):
    gene = global_session.query(Kinases).get(query)
    if gene is None:
        return redirect(url_for('user_message', query=query))
    gene_name = gene.Gene_symbol
    inhibitor_list = gene.inhibitors
    inh_number = len(inhibitor_list)
    return render_template("inhibitor_list.html", title='Inhibitor Search Results', inhibitor_list = inhibitor_list, gene_name=gene_name, inh_number=inh_number)

# The user it taken to this route when an inhibitor from the previous route is selected. 
# A new session is used. Background info is extracted for the inhibitor from the Inhibitor table.
@app.route("/search/inhibitors/summary/<path:inhibitor>")
def inhibitor_info(inhibitor):
    name = global_session.query(Inhibitors).get(inhibitor)
    url = name.Chemical_image
    formula = name.Chemical_Structure
    weight = name.Molecular_Weight
    kinase_list = name.kinases
    return render_template("inhibitor_info.html", title="Inhibitor", inhibitor=inhibitor, name=name, formula=formula, weight=weight, kinase_list=kinase_list, url=url)
    
# Create an extra route which displays an error message when no query matches are found.
@app.route("/error/<query>")
def user_message(query):
    return render_template("browse_message.html", title='No Matches', query=query)
    
if __name__ == '__main__':
    app.run(debug=False)
