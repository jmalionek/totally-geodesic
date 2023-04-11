from . import regina_util, surfaces, faces
import json
import regina
from sage.all import vector

reconstruction_example = {
 'name': 's783',
 'tri_used': 'gLLPQceeffefhuplalx',
 'vertex_surfaces': '{"N8":[0,4,0,0,0,0,0,0,0,0,0,3,0,0,4,1,0,0],"N3":[0,1,0,0,1,0,1,0,0,0,0,1,0,0,2,0,0,0],"N6":[1,0,0,2,0,0,0,1,0,0,0,0,1,0,0,1,0,0],"N7":[0,0,0,4,0,0,0,4,0,0,0,1,0,0,0,3,0,0]}',
 'vertex_genera': '{"N8":3,"N3":2,"N6":2,"N7":3}',
 'max_faces': '[{"dim":1,"verts":["N6","N7"]},{"dim":1,"verts":["N3","N8"]},{"dim":1,"verts":["N7","N8"]}]',
 'all_faces': '[{"dim":0,"verts":["N7"]},{"dim":0,"verts":["N3"]},{"dim":0,"verts":["N6"]},{"dim":0,"verts":["N8"]},{"dim":1,"verts":["N6","N7"]},{"dim":1,"verts":["N3","N8"]},{"dim":1,"verts":["N7","N8"]}]',
 'gen_func': '(-x^2 + 3*x)/(x^2 - 2*x + 1)',
}


def surfaces_to_zeroset_of_faces(surfaces):
    zeros = [faces.zero_set_standard(S) for S in surfaces]
    return frozenset.intersection(*zeros)

def reconstruct_faces(task):
    """
    >>> surfaces, LW = reconstruct_faces(reconstruction_example)
    """
    T = regina_util.as_regina(str(task['tri_used']))
    vertex_surfaces = dict()
    genera = json.loads(task['vertex_genera'])
    for label, quad_vector in json.loads(task['vertex_surfaces']).items():
        R = regina.NormalSurface(T, regina.NS_QUAD_CLOSED, quad_vector)
        assert R.isConnected() and R.eulerChar() == 2 - 2*genera[label] and label[0] == 'N'
        S = surfaces.NormalSurface(R, int(label[1:]))
        vertex_surfaces[label] = S

    all_faces = []
    for face in json.loads(task['all_faces']):
        its_surfaces = [vertex_surfaces[label] for label in face['verts']]
        zeros = surfaces_to_zeroset_of_faces(its_surfaces)
        all_faces.append(faces.AdmissibleFace(face['dim'], zeros, its_surfaces))

    zeros_of_max_faces = []
    max_faces = json.loads(task['max_faces'])
    for face in max_faces:
        its_surfaces = [vertex_surfaces[label] for label in face['verts']]
        zeros = surfaces_to_zeroset_of_faces(its_surfaces)
        zeros_of_max_faces.append(zeros)

    for face in all_faces:
        face.maximal = face.zero_set in zeros_of_max_faces

    assert len({face.zero_set for face in all_faces}) == len(all_faces)
    assert len([face for face in all_faces if face.maximal]) == len(max_faces)
    return vertex_surfaces, faces.FacesComplex(all_faces)

if __name__ == '__main__':
    import doctest
    doctest.testmod()

