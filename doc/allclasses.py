from sphinx.util.compat import Directive
from docutils import nodes
import pdb

class classlist(nodes.General, nodes.Element):
	pass

class ClasslistDirective(Directive):

	def run(self):
		return [classlist('')]

def process_classlist(app, doctree, fromdocname):
	env = app.builder.env

	para = nodes.paragraph()
	py = env.get_domain('py')
	classes = [_ for _ in py.get_objects() if _[2] == 'class']

	n = dict()
	for e in classes:
		newnode = nodes.reference('', '')
		innernode = nodes.Text(e[0])
		newnode['refdocname'] = e[3]
		newnode['refuri'] = app.builder.get_relative_uri(
                fromdocname, e[3])
		newnode['refuri'] += '#' + e[4]
		newnode.append(innernode)
		n[e[0].lower()] = newnode

	for key in sorted(n.iterkeys()):
		para += n[key]
		para += nodes.Text(' ')

	ctx = app.env.config['html_context']
	ctx['classlist'] = para
	for node in doctree.traverse(classlist):
		node.replace_self([para])
		continue

def add_classlist_handler(app):

	def _print_classlist(**kwargs):
		ctx = app.env.config['html_context']
		return app.builder.render_partial(ctx['classlist'])['fragment']

	ctx = app.env.config['html_context']
	if 'print_classlist' not in ctx:
		ctx['print_classlist'] = _print_classlist

def setup(app):
	app.add_node(classlist)
	app.add_directive('classlist', ClasslistDirective)
	app.connect('doctree-resolved', process_classlist)
	app.connect('builder-inited', add_classlist_handler)
