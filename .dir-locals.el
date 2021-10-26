((nil
  (eval . (let ((root (projectile-project-root)))
            (add-to-list
             (make-variable-buffer-local 'company-clang-arguments)
             (concat "-I" root "inst/include"))
            (add-to-list
             (make-variable-buffer-local 'flycheck-clang-include-path)
             (concat root "inst/include"))
            ))
  ))
