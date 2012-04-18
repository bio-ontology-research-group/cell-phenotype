import java.util.logging.Logger
import org.semanticweb.owlapi.apibinding.OWLManager
import org.semanticweb.owlapi.model.*
import org.semanticweb.owlapi.reasoner.*
import org.semanticweb.owlapi.profiles.*
import org.semanticweb.owlapi.util.*
import org.mindswap.pellet.KnowledgeBase
import org.mindswap.pellet.expressivity.*
import org.mindswap.pellet.*
import org.semanticweb.owlapi.io.*
import org.semanticweb.elk.owlapi.*
import org.semanticweb.owlapi.vocab.OWLRDFVocabulary

def ontfile = new File("gene_ontology_ext.obo")
def patofile = new File("quality.obo")
def chebifile = new File("chebi.obo")
def goxpfile = new File("go_xp_all.obo")

String formatClassNames(String s) {
  s=s.replace("<http://purl.obolibrary.org/obo/","")
  s=s.replace(">","")
  s=s.replace("_",":")
  s
}

def cli = new CliBuilder()
cli.with {
usage: 'Self'
  h longOpt:'help', 'this information'
  o longOpt:'output-file', 'output file',args:1, required:true
  //  t longOpt:'threads', 'number of threads', args:1
  //  k longOpt:'stepsize', 'steps before splitting jobs', arg:1
}

def opt = cli.parse(args)
if( !opt ) {
  //  cli.usage()
  return
}
if( opt.h ) {
  cli.usage()
  return
}

OWLOntologyManager manager = OWLManager.createOWLOntologyManager()

OWLDataFactory fac = manager.getOWLDataFactory()
OWLDataFactory factory = fac

def ontset = new TreeSet()
OWLOntology ont = manager.loadOntologyFromOntologyDocument(ontfile)
ontset.add(ont)
ont = manager.loadOntologyFromOntologyDocument(patofile)
ontset.add(ont)
ont = manager.loadOntologyFromOntologyDocument(chebifile)
ontset.add(ont)
ont = manager.createOntology(IRI.create("http://lc2.eu/temp.owl"), ontset)

OWLOntology outont = manager.createOntology(IRI.create("http://phenomebrowser.net/cellphenotype.owl"))
def onturi = "http://phenomebrowser.net/cellphenotype.owl#"

OWLReasonerFactory reasonerFactory = null

ConsoleProgressMonitor progressMonitor = new ConsoleProgressMonitor()
OWLReasonerConfiguration config = new SimpleConfiguration(progressMonitor)

OWLReasonerFactory f1 = new ElkReasonerFactory()
OWLReasoner reasoner = f1.createReasoner(ont,config)

OWLAnnotationProperty label = fac.getOWLAnnotationProperty(OWLRDFVocabulary.RDFS_LABEL.getIRI())

reasoner.precomputeInferences(InferenceType.CLASS_HIERARCHY)

def r = { String s ->
  factory.getOWLObjectProperty(IRI.create("http://phenomebrowser.net/ro2.owl#"+s))
}

def c = { String s ->
  factory.getOWLClass(IRI.create(onturi+s))
}

def id2class = [:] // maps a name to an OWLClass
ont.getClassesInSignature(true).each {
  def aa = it.toString()
  aa = formatClassNames(aa)
  if (id2class[aa] != null) {
  } else {
    id2class[aa] = it
  }
}


def addAnno = {resource, prop, cont ->
  OWLAnnotation anno = factory.getOWLAnnotation(
    factory.getOWLAnnotationProperty(prop.getIRI()),
    factory.getOWLTypedLiteral(cont))
  def axiom = factory.getOWLAnnotationAssertionAxiom(resource.getIRI(),
                                                     anno)
  manager.addAxiom(outont,axiom)
}


/* Parse the GO-XP file first to create map between CC and BP */
def ccs = reasoner.getSubClasses(id2class["GO:0005575"], false).getFlattened()
def bps = reasoner.getSubClasses(id2class["GO:0009987"], false).getFlattened()
def regs = reasoner.getSubClasses(id2class["GO:0050789"], false).getFlattened()

def inputmap = [:]
def outputmap = [:]

def oboid = ""
def bp2cc = [:]
def cc2bp = [:]
goxpfile.eachLine { line ->
  if (line.startsWith("id:")) {
    oboid = line.substring(4).trim()
  }
  if (line.startsWith("intersection_of: OBO_REL:") && line.indexOf("GO:")>-1) {
    def sid = ""
    if (line.indexOf("!")>-1) {
      sid = line.substring(0, line.indexOf("!")).trim()
    } else {
      sid = line.trim()
    }
    sid = sid.substring(sid.indexOf(" ")).trim()
    sid = sid.substring(sid.indexOf(" ")).trim()
    if ((id2class[oboid] in bps) && (id2class[sid] in ccs)) {
      if (bp2cc[oboid]==null) {
	bp2cc[oboid] = new TreeSet()
      }
      bp2cc[oboid].add(sid)
      if (cc2bp[sid] == null) {
	cc2bp[sid] = new TreeSet()
      }
      cc2bp[sid].add(oboid)
    }
  }
  if (line.startsWith("intersection_of: OBO_REL:has_input")) {
    def sid = ""
    if (line.indexOf("!")>-1) {
      sid = line.substring(0, line.indexOf("!")).trim()
    } else {
      sid = line.trim()
    }
    sid = sid.substring(sid.indexOf(" ")).trim()
    sid = sid.substring(sid.indexOf(" ")).trim()
    if (inputmap[oboid] == null) {
      inputmap[oboid] = new TreeSet()
    }
    inputmap[oboid].add(sid)
  }  
  if (line.startsWith("intersection_of: OBO_REL:has_output")) {
    def sid = ""
    if (line.indexOf("!")>-1) {
      sid = line.substring(0, line.indexOf("!")).trim()
    } else {
      sid = line.trim()
    }
    sid = sid.substring(sid.indexOf(" ")).trim()
    sid = sid.substring(sid.indexOf(" ")).trim()
    if (outputmap[oboid] == null) {
      outputmap[oboid] = new TreeSet()
    }
    outputmap[oboid].add(sid)
  }  
}

addAnno(c("C3PO:000000000"),OWLRDFVocabulary.RDFS_LABEL,"cellular component phenotype")
addAnno(c("C3PO:999999999"),OWLRDFVocabulary.RDFS_LABEL,"cellular process phenotype")

/* Start with Cell Components */
ccs.each { cc ->
  def id = formatClassNames(cc.toString()).substring(3) // id of GO class
  if (! id.startsWith(":")) {
    def nid = "C3PO:00$id" // create new C3PO class with prefix "00"
    def cl = c(nid)
    manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl, c("C3PO:000000000")))
    def name = null // label of GO Class
    cc.getAnnotations(ont, label).each { name = it.getValue().getLiteral() }
    addAnno(cl,OWLRDFVocabulary.RDFS_LABEL,"$name phenotype")
    addAnno(cl,OWLRDFVocabulary.RDF_DESCRIPTION,"Any observable characteristic of $name.")
    manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
		       cl,
		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
						      fac.getOWLObjectSomeValuesFrom(
							r("has-part"), fac.getOWLObjectSomeValuesFrom(
							  r("part-of"), cc)))))
  
    // def cl2 = c("C3PO:01$id") 
    // manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl2, cl))
    // addAnno(cl2,OWLRDFVocabulary.RDFS_LABEL,"Normal $name phenotype")
    // addAnno(cl2,OWLRDFVocabulary.RDF_DESCRIPTION,"The observable characteristics of $name are normal.")

    // def cl3 = c("C3PO:02$id") 
    // manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl3, cl))
    // addAnno(cl3,OWLRDFVocabulary.RDFS_LABEL,"Abnormal $name phenotype")
    // addAnno(cl3,OWLRDFVocabulary.RDF_DESCRIPTION,"Some observable characteristics of $name are abnormal.")
    // manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 		       cl3,
    // 		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 						      fac.getOWLObjectSomeValuesFrom(
    // 							r("has-part"), fac.getOWLObjectIntersectionOf(
    // 							  fac.getOWLObjectSomeValuesFrom(
    // 							    r("part-of"), cc),
    // 							  fac.getOWLObjectSomeValuesFrom(
    // 							    r("has-quality"), id2class["PATO:0000001"]))))))

    // def cl4 = c("C3PO:03$id") 
    // manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl4, cl3))
    // addAnno(cl4,OWLRDFVocabulary.RDFS_LABEL,"Absence of $name")
    // addAnno(cl4,OWLRDFVocabulary.RDF_DESCRIPTION,"The complete absence of $name. No $name is present in the organism.")

    def cl5 = c("C3PO:54$id") 
    manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl5, cl))
    addAnno(cl5,OWLRDFVocabulary.RDFS_LABEL,"$name morphology phenotype")
    addAnno(cl5,OWLRDFVocabulary.RDF_DESCRIPTION,"Morphology phenotype of $name.")
    manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
		       cl5,
		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
						      fac.getOWLObjectSomeValuesFrom(
							r("has-part"), fac.getOWLObjectIntersectionOf(
							  fac.getOWLObjectSomeValuesFrom(r("part-of"), cc), 
							  fac.getOWLObjectSomeValuesFrom(
							    r("has-quality"), id2class["PATO:0000051"]))))))

    def cl6 = c("C3PO:55$id") 
    manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl6, cl))
    addAnno(cl6,OWLRDFVocabulary.RDFS_LABEL,"$name physiology phenotype")
    addAnno(cl6,OWLRDFVocabulary.RDF_DESCRIPTION,"Physiology phenotype of $name. The physiology of $name includes $name's functions and functionings.")
    manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
		       cl6,
		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
						      fac.getOWLObjectSomeValuesFrom(
							r("has-part"), fac.getOWLObjectIntersectionOf(
							  cc, fac.getOWLObjectSomeValuesFrom(
							    r("has-quality"), id2class["PATO:0001509"]))))))
  }
}

bps.each { bp ->
  def id = formatClassNames(bp.toString()).substring(3) // id of GO class
  if (! id.startsWith(":")) {
    def nid = "C3PO:10$id" // create new C3PO class with prefix "00"
    def cl = c(nid)
    manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl, c("C3PO:999999999")))
    def name = null // label of GO Class
    bp.getAnnotations(ont, label).each { name = it.getValue().getLiteral() }
    addAnno(cl,OWLRDFVocabulary.RDFS_LABEL,"$name phenotype")
    addAnno(cl,OWLRDFVocabulary.RDF_DESCRIPTION,"An observable characteristic of processes of the type $name.")

    manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
		       cl,
		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
						      fac.getOWLObjectSomeValuesFrom(
							r("has-part"), fac.getOWLObjectSomeValuesFrom(
							  r("participates-in"), 
							  fac.getOWLObjectSomeValuesFrom(
							    r("part-of"), bp))))))

  
    // def cl2 = c("C3PO:11$id") 
    // manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl2, cl))
    // addAnno(cl2,OWLRDFVocabulary.RDFS_LABEL,"Normal $name phenotype")
    // addAnno(cl2,OWLRDFVocabulary.RDF_DESCRIPTION,"The observable characteristics of $name processes are normal.")

    // def cl3 = c("C3PO:12$id") 
    // manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl3, cl))
    // addAnno(cl3,OWLRDFVocabulary.RDFS_LABEL,"Abnormal $name phenotype")
    // addAnno(cl3,OWLRDFVocabulary.RDF_DESCRIPTION,"The observable characteristics of $name processes are abnormal.")
    /* This does not work, because we need to use a union to say that it is either the process that is abnormal of the regulation of the process */
    // manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 		     cl3,
    // 		     fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 						    fac.getOWLObjectSomeValuesFrom(
    // 						      r("has-part"), 
    // 						      fac.getOWLObjectSomeValuesFrom(
    // 							r("participates-in"), 
    // 							fac.getOWLObjectIntersectionOf(
    // 							  fac.getOWLObjectSomeValuesFrom(
    // 							    r("part-of"), bp),
    // 							  fac.getOWLObjectSomeValuesFrom(
    // 							    r("has-quality"), id2class["PATO:0000001"])))))))

    if (bp2cc["GO:"+id]!=null) {
      bp2cc["GO:"+id].each { cid2 ->
	cid2 = cid2.substring(3)
	def cccl = c("C3PO:55$cid2")
	manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl,cccl))
      }
    }

    // def cl4 = c("C3PO:13$id") 
    // manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl4, cl))
    // addAnno(cl4,OWLRDFVocabulary.RDFS_LABEL,"Absence of $name")
    // addAnno(cl4,OWLRDFVocabulary.RDF_DESCRIPTION,"The complete absence of $name processes.")

    def cl5 = c("C3PO:14$id") 
    manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl5, cl))
    addAnno(cl5,OWLRDFVocabulary.RDFS_LABEL,"phenotype of single occurrence of $name")
    addAnno(cl5,OWLRDFVocabulary.RDF_DESCRIPTION,"Phenotype of single $name occurrences. Examples include (extended) durations, (abnormal) course of the process, or abnormalities of process participants.")
    manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
		       cl5,
		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
						      fac.getOWLObjectSomeValuesFrom(
							r("has-part"), 
							fac.getOWLObjectSomeValuesFrom(
							  r("participates-in"), 
							  fac.getOWLObjectIntersectionOf(
							    fac.getOWLObjectSomeValuesFrom(
							      r("part-of"), bp),
							    fac.getOWLObjectSomeValuesFrom(
							      r("has-quality"), id2class["PATO:0000001"])))))))

    def cla = c("C3PO:55$id") 
    manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cla, cl5))
    addAnno(cla,OWLRDFVocabulary.RDFS_LABEL,"duration phenotype of $name")
    addAnno(cla,OWLRDFVocabulary.RDF_DESCRIPTION,"The duration of $name. Duration can have a value, or be increased or decreased.")
    manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    		       cla,
    		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    						      fac.getOWLObjectSomeValuesFrom(
    							r("has-part"), 
    							fac.getOWLObjectSomeValuesFrom(
    							  r("participates-in"), 
    							  fac.getOWLObjectIntersectionOf(
    							    bp,
    							    fac.getOWLObjectSomeValuesFrom(
    							      r("has-quality"), id2class["PATO:0001309"])))))))
  
    // cla = c("C3PO:56$id") 
    // manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cla, cl3))
    // addAnno(cla,OWLRDFVocabulary.RDFS_LABEL,"Decreased duration of $name")
    // addAnno(cla,OWLRDFVocabulary.RDF_DESCRIPTION,"The duration of $name is decreased compared to normal.")
    // manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 		       cla,
    // 		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 						      fac.getOWLObjectSomeValuesFrom(
    // 							r("has-part"), 
    // 							fac.getOWLObjectSomeValuesFrom(
    // 							  r("participates-in"), 
    // 							  fac.getOWLObjectIntersectionOf(
    // 							    bp,
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-quality"), id2class["PATO:0000499"])))))))
  
    def cl6 = c("C3PO:15$id") 
    manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl6, cl))
    addAnno(cl6,OWLRDFVocabulary.RDFS_LABEL,"phenotype of $name regulation")
    addAnno(cl6,OWLRDFVocabulary.RDF_DESCRIPTION,"Phenotype of processes which regulate occurrences of $name. Examples include (increased) frequency of occurrence, (late, early) onset, or abnormalities of process participants throughout the total course of the regulation process.")
    manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
		       cl6,
		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
						      fac.getOWLObjectSomeValuesFrom(
							r("has-part"),
							fac.getOWLObjectSomeValuesFrom(
							  r("participates-in"), 
							  fac.getOWLObjectIntersectionOf(
							    fac.getOWLObjectSomeValuesFrom(
							      r("regulates"), bp),
							    fac.getOWLObjectSomeValuesFrom(
							      r("has-quality"), 
							      fac.getOWLObjectIntersectionOf(
								id2class["PATO:0000001"], fac.getOWLObjectSomeValuesFrom(r("towards"), bp)))))))))
		   
    def cl7 = c("C3PO:16$id") 
    manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl7, cl6))
    addAnno(cl7,OWLRDFVocabulary.RDFS_LABEL,"frequency of $name occurrences")
    addAnno(cl7,OWLRDFVocabulary.RDF_DESCRIPTION,"The frequency of occurrences of $name.")
    manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    		       cl7,
    		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    						      fac.getOWLObjectSomeValuesFrom(
    							r("has-part"), 
    							fac.getOWLObjectSomeValuesFrom(
    							  r("participates-in"), 
    							  fac.getOWLObjectIntersectionOf(
    							    fac.getOWLObjectSomeValuesFrom(
    							      r("regulates"), bp),
    							    fac.getOWLObjectSomeValuesFrom(
    							      r("has-quality"), 
    							      fac.getOWLObjectIntersectionOf(
    								id2class["PATO:0000044"], fac.getOWLObjectSomeValuesFrom(r("towards"), bp)))))))))

    // def cl8 = c("C3PO:17$id") 
    // manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl8, cl6))
    // addAnno(cl8,OWLRDFVocabulary.RDFS_LABEL,"Decreased frequency of $name")
    // addAnno(cl8,OWLRDFVocabulary.RDF_DESCRIPTION,"A regulatory abnormality in which the frequency of occurrences of $name is decreased.")
    // manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 		       cl8,
    // 		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 						      fac.getOWLObjectSomeValuesFrom(
    // 							r("has-part"), 
    // 							fac.getOWLObjectSomeValuesFrom(
    // 							  r("participates-in"), 
    // 							  fac.getOWLObjectIntersectionOf(
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("regulates"), bp),
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-quality"), 
    // 							      fac.getOWLObjectIntersectionOf(
    // 								id2class["PATO:0000381"], fac.getOWLObjectSomeValuesFrom(r("towards"), bp)))))))))

    def clx = c("C3PO:90$id") 
    manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(clx, cl6))
    addAnno(clx,OWLRDFVocabulary.RDFS_LABEL,"onset phenotype of $name")
    addAnno(clx,OWLRDFVocabulary.RDF_DESCRIPTION,"The onset of $name.")
    manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    		       clx,
    		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    						      fac.getOWLObjectSomeValuesFrom(
    							r("has-part"), 
    							fac.getOWLObjectSomeValuesFrom(
    							  r("participates-in"), 
    							  fac.getOWLObjectIntersectionOf(
    							    fac.getOWLObjectSomeValuesFrom(
    							      r("regulates"), bp),
    							    fac.getOWLObjectSomeValuesFrom(
    							      r("has-quality"), 
    							      fac.getOWLObjectIntersectionOf(
    								id2class["PATO:0002325"], fac.getOWLObjectSomeValuesFrom(r("towards"), bp)))))))))

    // def cly = c("C3PO:91$id") 
    // manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cly, cl6))
    // addAnno(cly,OWLRDFVocabulary.RDFS_LABEL,"Early onset of $name")
    // addAnno(cly,OWLRDFVocabulary.RDF_DESCRIPTION,"A regulatory abnormality in which $name starts earlier than normal.")
    // manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 		       cly,
    // 		       fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 						      fac.getOWLObjectSomeValuesFrom(
    // 							r("has-part"), 
    // 							fac.getOWLObjectSomeValuesFrom(
    // 							  r("participates-in"), 
    // 							  fac.getOWLObjectIntersectionOf(
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("regulates"), bp),
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-quality"), 
    // 							      fac.getOWLObjectIntersectionOf(
    // 								id2class["PATO:0000694"], fac.getOWLObjectSomeValuesFrom(r("towards"), bp)))))))))

    // if (outputmap["GO:$id"]!=null) {
    //   def output = outputmap["GO:$id"]
    //   def oname = null
    //   output.each { outp ->
    // 	id2class[outp]?.getAnnotations(ont, label).each { oname = it.getValue().getLiteral() }
    // 	if (oname!=null) {
    // 	  def cl9 = c("C3PO:18$id")
    // 	  manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl9, cl6))
    // 	  addAnno(cl9,OWLRDFVocabulary.RDFS_LABEL,"Decreased mass of $oname as output in regulation of $name")
    // 	  addAnno(cl9,OWLRDFVocabulary.RDF_DESCRIPTION,"The total mass of $oname that is used as output throughout one or more $name processes is decreased.")
    // 	  manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 			     cl9,
    // 			     fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-part"), 
    // 							      fac.getOWLObjectSomeValuesFrom(
    // 								r("participates-in"),
    // 								fac.getOWLObjectIntersectionOf(
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("regulates"), bp),
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("has-output"),
    // 								    fac.getOWLObjectIntersectionOf(
    // 								      id2class[outp],
    // 								      fac.getOWLObjectSomeValuesFrom(
    // 									r("has-quality"), 
    // 									id2class["PATO:0001562"])))))))))
    // 	  def cl10 = c("C3PO:19$id") 
    // 	  manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl10, cl6))
    // 	  addAnno(cl10,OWLRDFVocabulary.RDFS_LABEL,"Increased mass of $oname as output in regulation of $name")
    // 	  addAnno(cl10,OWLRDFVocabulary.RDF_DESCRIPTION,"The total mass of $oname that is used as output throughout one or more $name processes is increased.")
    // 	  manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 			     cl10,
    // 			     fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-part"), 
    // 							      fac.getOWLObjectSomeValuesFrom(
    // 								r("participates-in"),
    // 								fac.getOWLObjectIntersectionOf(
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("regulates"), bp),
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("has-output"),
    // 								    fac.getOWLObjectIntersectionOf(
    // 								      id2class[outp],
    // 								      fac.getOWLObjectSomeValuesFrom(
    // 									r("has-quality"), 
    // 									id2class["PATO:0001563"])))))))))
    // 	  def cl11 = c("C3PO:20$id") 
    // 	  manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl11, cl5))
    // 	  addAnno(cl11,OWLRDFVocabulary.RDFS_LABEL,"Increased mass of $oname as output in $name")
    // 	  addAnno(cl11,OWLRDFVocabulary.RDF_DESCRIPTION,"The mass of $oname that is used as output in a single $name is increased.")
    // 	  manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 			     cl11,
    // 			     fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-part"), 
    // 							      fac.getOWLObjectSomeValuesFrom(
    // 								r("participates-in"),
    // 								fac.getOWLObjectIntersectionOf(
    // 								  bp,
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("has-output"),
    // 								    fac.getOWLObjectIntersectionOf(
    // 								      id2class[outp],
    // 								      fac.getOWLObjectSomeValuesFrom(
    // 									r("has-quality"), 
    // 									id2class["PATO:0001563"])))))))))
    // 	  def cl12 = c("C3PO:21$id") 
    // 	  manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl12, cl5))
    // 	  addAnno(cl12,OWLRDFVocabulary.RDFS_LABEL,"Decreased mass of $oname as output in $name")
    // 	  addAnno(cl12,OWLRDFVocabulary.RDF_DESCRIPTION,"The mass of $oname that is used as output in a single $name is decreased.")
    // 	  manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 			     cl12,
    // 			     fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-part"), 
    // 							      fac.getOWLObjectSomeValuesFrom(
    // 								r("participates-in"),
    // 								fac.getOWLObjectIntersectionOf(
    // 								  bp,
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("has-output"),
    // 								    fac.getOWLObjectIntersectionOf(
    // 								      id2class[outp],
    // 								      fac.getOWLObjectSomeValuesFrom(
    // 									r("has-quality"), 
    // 									id2class["PATO:0001562"])))))))))
    // 	}
    //   }
    // }

    // if (inputmap["GO:$id"]!=null) {
    //   def output = inputmap["GO:$id"]
    //   def oname = null
    //   output.each { outp ->
    // 	id2class[outp]?.getAnnotations(ont, label).each { oname = it.getValue().getLiteral() }
    // 	if (oname!=null) {
    // 	  def cl9 = c("C3PO:22$id")
    // 	  manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl9, cl6))
    // 	  addAnno(cl9,OWLRDFVocabulary.RDFS_LABEL,"Decreased mass of $oname as input in regulation of $name")
    // 	  addAnno(cl9,OWLRDFVocabulary.RDF_DESCRIPTION,"The total mass of $oname that is used as input throughout one or more $name processes is decreased.")
    // 	  manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 			     cl9,
    // 			     fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-part"), 
    // 							      fac.getOWLObjectSomeValuesFrom(
    // 								r("participates-in"),
    // 								fac.getOWLObjectIntersectionOf(
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("regulates"), bp),
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("has-input"),
    // 								    fac.getOWLObjectIntersectionOf(
    // 								      id2class[outp],
    // 								      fac.getOWLObjectSomeValuesFrom(
    // 									r("has-quality"), 
    // 									id2class["PATO:0001562"])))))))))
    // 	  def cl10 = c("C3PO:23$id") 
    // 	  manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl10, cl6))
    // 	  addAnno(cl10,OWLRDFVocabulary.RDFS_LABEL,"Increased mass of $oname as input in regulation of $name")
    // 	  addAnno(cl10,OWLRDFVocabulary.RDF_DESCRIPTION,"The total mass of $oname that is used as input throughout one or more $name processes is increased.")
    // 	  manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 			     cl10,
    // 			     fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-part"), 
    // 							      fac.getOWLObjectSomeValuesFrom(
    // 								r("participates-in"),
    // 								fac.getOWLObjectIntersectionOf(
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("regulates"), bp),
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("has-input"),
    // 								    fac.getOWLObjectIntersectionOf(
    // 								      id2class[outp],
    // 								      fac.getOWLObjectSomeValuesFrom(
    // 									r("has-quality"), 
    // 									id2class["PATO:0001563"])))))))))
    // 	  def cl11 = c("C3PO:24$id") 
    // 	  manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl11, cl5))
    // 	  addAnno(cl11,OWLRDFVocabulary.RDFS_LABEL,"Increased mass of $oname as input in $name")
    // 	  addAnno(cl11,OWLRDFVocabulary.RDF_DESCRIPTION,"The mass of $oname that is used as input in a single $name is increased.")
    // 	  manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 			     cl11,
    // 			     fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-part"), 
    // 							      fac.getOWLObjectSomeValuesFrom(
    // 								r("participates-in"),
    // 								fac.getOWLObjectIntersectionOf(
    // 								  bp,
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("has-input"),
    // 								    fac.getOWLObjectIntersectionOf(
    // 								      id2class[outp],
    // 								      fac.getOWLObjectSomeValuesFrom(
    // 									r("has-quality"), 
    // 									id2class["PATO:0001563"])))))))))
    // 	  def cl12 = c("C3PO:25$id") 
    // 	  manager.addAxiom(outont, factory.getOWLSubClassOfAxiom(cl12, cl5))
    // 	  addAnno(cl12,OWLRDFVocabulary.RDFS_LABEL,"Decreased mass of $oname as input in $name")
    // 	  addAnno(cl12,OWLRDFVocabulary.RDF_DESCRIPTION,"The mass of $oname that is used as input in a single $name is decreased.")
    // 	  manager.addAxiom(outont, fac.getOWLEquivalentClassesAxiom(
    // 			     cl12,
    // 			     fac.getOWLObjectSomeValuesFrom(r("phenotype-of"),
    // 							    fac.getOWLObjectSomeValuesFrom(
    // 							      r("has-part"), 
    // 							      fac.getOWLObjectSomeValuesFrom(
    // 								r("participates-in"),
    // 								fac.getOWLObjectIntersectionOf(
    // 								  bp,
    // 								  fac.getOWLObjectSomeValuesFrom(
    // 								    r("has-input"),
    // 								    fac.getOWLObjectIntersectionOf(
    // 								      id2class[outp],
    // 								      fac.getOWLObjectSomeValuesFrom(
    // 									r("has-quality"), 
    // 									id2class["PATO:0001562"])))))))))
    // 	}

  }
}

manager.addAxiom(
  outont, fac.getOWLEquivalentObjectPropertiesAxiom(
    r("has-part"), fac.getOWLObjectProperty(IRI.create("http://purl.obolibrary.org/obo/GENE_ONTOLOGY_has_part"))))

manager.addAxiom(
  outont, fac.getOWLEquivalentObjectPropertiesAxiom(
    r("part-of"), fac.getOWLObjectProperty(IRI.create("http://purl.obolibrary.org/obo/GENE_ONTOLOGY_part_of"))))

manager.addAxiom(
  outont, fac.getOWLEquivalentObjectPropertiesAxiom(
    r("part-of"), fac.getOWLObjectProperty(IRI.create("http://purl.obolibrary.org/obo/BFO_0000050"))))

manager.addAxiom(
  outont, fac.getOWLEquivalentObjectPropertiesAxiom(
    r("regulates"), fac.getOWLObjectProperty(IRI.create("http://purl.obolibrary.org/obo/GENE_ONTOLOGY_regulates"))))


manager.addAxiom(outont, fac.getOWLTransitiveObjectPropertyAxiom(r("has-part")))
manager.addAxiom(outont, fac.getOWLTransitiveObjectPropertyAxiom(r("part-of")))
manager.addAxiom(outont, fac.getOWLTransitiveObjectPropertyAxiom(r("regulates")))
manager.addAxiom(outont, fac.getOWLReflexiveObjectPropertyAxiom(r("has-part")))
manager.addAxiom(outont, fac.getOWLReflexiveObjectPropertyAxiom(r("part-of")))
manager.addAxiom(outont, fac.getOWLReflexiveObjectPropertyAxiom(r("regulates")))
OWLImportsDeclaration importDecl1 = fac.getOWLImportsDeclaration(IRI.create("http://purl.obolibrary.org/obo/go.owl"))
manager.applyChange(new AddImport(outont, importDecl1))
importDecl1 = fac.getOWLImportsDeclaration(IRI.create("http://purl.obolibrary.org/obo/pato.owl"))
manager.applyChange(new AddImport(outont, importDecl1))
importDecl1 = fac.getOWLImportsDeclaration(IRI.create("http://purl.obolibrary.org/obo/chebi.owl"))
manager.applyChange(new AddImport(outont, importDecl1))

manager.saveOntology(outont, IRI.create("file:"+opt.o))
System.exit(0)
